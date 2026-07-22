"""Batched local amplitude scaling.
"""
import numpy as np
import torch


def round_up_proper(x):
    """scaling_tools' round_up_proper: round, nudged so exact .5 rounds up."""
    return int(np.round(x + 1e-5))


def calculate_frequency_map(shape):
    """Per-voxel frequency magnitude over the rfft grid (FDRutil.calculate_frequency_map)."""
    fi = np.fft.fftfreq(shape[0], 1.0)
    fj = np.fft.fftfreq(shape[1], 1.0)
    fk = np.fft.rfftfreq(shape[2], 1.0)
    grid = (fi.size, fj.size, fk.size)
    fmi = (fi * fi)[:, None, None] * np.ones(grid)
    fmj = (fj * fj)[None, :, None] * np.ones(grid)
    fmk = (fk * fk)[None, None, :] * np.ones(grid)
    return np.sqrt(fmi + fmj + fmk)


def pick_device(device=None):
    if device is not None:
        return torch.device(device)
    if torch.cuda.is_available():
        return torch.device("cuda")
    if torch.backends.mps.is_available():
        return torch.device("mps")
    return torch.device("cpu")


class WindowedScaler:
    """Scales a batch of target windows so each one's radial amplitude profile matches its
    reference window, then returns the central voxel of the (inverse-FFT'd) scaled window.

    Everything shape/frequency related is precomputed once, in __init__.
    """

    def __init__(self, wn, device="cpu", dtype=torch.float32):
        self.wn = wn
        self.device = device
        self.dtype = dtype
        self.cdtype = torch.complex128 if dtype == torch.float64 else torch.complex64

        self.fshape = (wn, wn, wn // 2 + 1)                 # rfftn output shape
        self.central_pix = round_up_proper(wn / 2.0)        # 13 for wn=25

        # shell index of every rfft voxel -- identical for every window, so computed once
        x, y, z = np.indices(self.fshape)
        radii = np.sqrt(x ** 2 + y ** 2 + z ** 2).astype(int)
        self.n_prof = wn // 2 + 1                            # profile truncated to [0 : wn/2+1]
        self.n_shell = int(radii.max()) + 1
        self.radii = torch.tensor(radii.ravel(), device=device)
        counts = np.bincount(radii.ravel(), minlength=self.n_shell).astype(np.float64)
        self.shell_counts = torch.tensor(counts, device=device, dtype=dtype)

        # interpolation axes for applying the scale factor
        self.frequencies = torch.tensor(np.fft.rfftfreq(wn), device=device, dtype=dtype)
        fmap = calculate_frequency_map((wn, wn, wn))
        self.frequency_map = torch.tensor(fmap.ravel(), device=device, dtype=dtype)

    def _radial_profile(self, F_abs):
        """(B, Nvox) amplitudes -> (B, n_prof) mean amplitude per shell."""
        B = F_abs.shape[0]
        prof = torch.zeros(B, self.n_shell, device=self.device, dtype=self.dtype)
        prof.index_add_(1, self.radii, F_abs)               # sum each voxel into its shell
        prof = prof / self.shell_counts                     # bincount average
        return prof[:, :self.n_prof]

    def _interp(self, scale_factors):
        """Batched np.interp(frequency_map, frequencies, scale_factors), ends clamped."""
        xp, fp, x = self.frequencies, scale_factors, self.frequency_map
        idx = torch.searchsorted(xp, x).clamp(1, xp.numel() - 1)
        x0, x1 = xp[idx - 1], xp[idx]
        w = ((x - x0) / (x1 - x0)).clamp(0, 1)
        f0, f1 = fp[:, idx - 1], fp[:, idx]
        return f0 + (f1 - f0) * w                           # (B, Nvox)

    def _scale_maps(self, reference_windows, target_windows):
        """FFT -> profile -> scale factor -> apply -> inverse FFT.

        Returns the whole scaled target cubes, (B, wn, wn, wn) real. The batch axis is left
        untouched by every transform; each window is an independent problem.
        """
        tar = torch.as_tensor(target_windows, device=self.device, dtype=self.dtype)
        ref = torch.as_tensor(reference_windows, device=self.device, dtype=self.dtype)
        B = tar.shape[0]

        Ftar = torch.fft.rfftn(tar, dim=(1, 2, 3), norm="ortho")   # (B, wn, wn, wn//2+1)
        Fref = torch.fft.rfftn(ref, dim=(1, 2, 3), norm="ortho")

        tar_prof = self._radial_profile(Ftar.abs().reshape(B, -1))
        ref_prof = self._radial_profile(Fref.abs().reshape(B, -1))

        scale = ref_prof / tar_prof                          # |ref| / |target| per shell
        scale = torch.where(torch.isfinite(scale), scale, torch.zeros_like(scale))

        scaling_map = self._interp(scale).reshape(B, *self.fshape).to(self.cdtype)
        scaled_fft = scaling_map * Ftar
        return torch.fft.irfftn(scaled_fft, s=(self.wn,) * 3, dim=(1, 2, 3), norm="ortho")

    def scale_central(self, reference_windows, target_windows):
        """(B, wn, wn, wn) x2 -> (B,) central voxel of each scaled target window."""
        scaled = self._scale_maps(reference_windows, target_windows)
        c = self.central_pix
        return scaled[:, c, c, c]


def _gather_windows(volume, corners, wn):
    """(N, wn, wn, wn) stack of windows, one per corner (k, j, i). Preserves volume dtype."""
    out = np.empty((len(corners), wn, wn, wn), dtype=volume.dtype)
    for n, (k, j, i) in enumerate(corners):
        out[n] = volume[k:k + wn, j:j + wn, i:i + wn]
    return out


def local_amplitude_scaling(reference_map, target_map, corner_positions, window_size,
                            device=None, chunk=4096, dtype=torch.float32):
    """Scale the local amplitudes of target_map to match reference_map, at each masked voxel.

    Returns `sharpened_vals`: the central scaled voxel for each entry, in the input order, so
    it lines up with `masked_indices` for put_scaled_voxels_back_in_original_volume.
    """
    from tqdm import tqdm
    
    np_dtype = np.float64 if dtype == torch.float64 else np.float32
    reference_map = np.asarray(reference_map, dtype=np_dtype)
    target_map = np.asarray(target_map, dtype=np_dtype)
    wn = int(window_size)
    centres = np.asarray(corner_positions)

    # centre -> corner, exactly as the old loop: round_up_proper(centre - wn/2)
    corners = np.round(centres - wn / 2.0 + 1e-5).astype(int)
    corners = np.clip(corners, 0, np.array(target_map.shape) - wn)   # keep the window in-bounds

    device = pick_device(device)
    print(f"Device is: {device}")
    scaler = WindowedScaler(wn, device=str(device), dtype=dtype)

    sharpened = np.empty(len(centres), dtype=np_dtype)
    for start in tqdm(range(0, len(centres), chunk), desc="LocScale"):
        sel = corners[start:start + chunk]
        tw = _gather_windows(target_map, sel, wn)
        rw = _gather_windows(reference_map, sel, wn)
        vals = scaler.scale_central(rw, tw)
        sharpened[start:start + len(sel)] = vals.cpu().numpy()

    return sharpened
