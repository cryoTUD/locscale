"""GPU-accelerated windowed amplitude scaling.
"""
import numpy as np
import torch


def central_pix(wn):
    """LocScale's round_up_proper(wn/2.0). 13 for wn=25, i.e. one past the geometric centre."""
    return int(np.round(wn / 2.0 + 1e-5))


def calculate_frequency_map(map):
    """Per-voxel frequency magnitude over the rfft grid (FDRutil.calculate_frequency_map)."""
    s = map.shape
    fi = np.fft.fftfreq(s[0], 1.0)
    fj = np.fft.fftfreq(s[1], 1.0)
    fk = np.fft.rfftfreq(s[2], 1.0)
    shape = (fi.size, fj.size, fk.size)
    fmi = (fi * fi)[:, None, None] * np.ones(shape)
    fmj = (fj * fj)[None, :, None] * np.ones(shape)
    fmk = (fk * fk)[None, None, :] * np.ones(shape)
    return np.sqrt(fmi + fmj + fmk)


class WindowedScaler:
    """Batched windowed scaling. Everything shape/frequency related is precomputed once."""

    def __init__(self, wn, device="cpu", dtype=torch.float32):
        self.wn = wn
        self.device = device
        self.dtype = dtype
        self.cdtype = torch.complex128 if dtype == torch.float64 else torch.complex64

        self.fshape = (wn, wn, wn // 2 + 1)
        self.central_pix = central_pix(wn)

        x, y, z = np.indices(self.fshape)
        radii = np.sqrt(x**2 + y**2 + z**2).astype(int)
        self.n_prof = wn // 2 + 1
        self.n_shell = int(radii.max()) + 1
        self.radii = torch.tensor(radii.ravel(), device=device)
        counts = np.bincount(radii.ravel(), minlength=self.n_shell).astype(np.float64)
        self.shell_counts = torch.tensor(counts, device=device, dtype=dtype)

        frequencies = torch.tensor(np.fft.rfftfreq(wn), device=device, dtype=dtype)
        frequency_map = torch.tensor(calculate_frequency_map(np.zeros((wn, wn, wn))).ravel(),
                                     device=device, dtype=dtype)

        # centre-voxel weight, so no inverse FFT is needed
        c = self.central_pix
        delta = np.zeros((wn, wn, wn)); delta[c, c, c] = 1.0
        Bd = np.fft.rfftn(delta, norm="ortho")
        mult = np.full(self.fshape, 2.0); mult[..., 0] = 1.0
        if wn % 2 == 0:
            mult[..., -1] = 1.0
        center_weight = torch.tensor(np.conj(Bd) * mult, device=device, dtype=self.cdtype)

        # fold np.interp(frequency_map, frequencies, scale) into shell-space scatter targets
        idx = torch.searchsorted(frequencies, frequency_map).clamp(1, frequencies.numel() - 1)
        f0, f1 = frequencies[idx - 1], frequencies[idx]
        w = ((frequency_map - f0) / (f1 - f0)).clamp(0, 1)
        self._idx0 = (idx - 1).contiguous()
        self._idx1 = idx.contiguous()
        Wflat = center_weight.reshape(-1)
        self._Wa = (Wflat * (1 - w).to(self.cdtype)).contiguous()
        self._Wb = (Wflat * w.to(self.cdtype)).contiguous()

    def _radial_profile(self, F_abs):
        """(B, Nvox) amplitudes -> (B, n_prof) mean amplitude per shell."""
        B = F_abs.shape[0]
        prof = torch.zeros(B, self.n_shell, device=self.device, dtype=self.dtype)
        prof.index_add_(1, self.radii, F_abs)
        prof = prof / self.shell_counts
        return prof[:, :self.n_prof]

    def scale_windows(self, target_windows, reference_windows):
        """(B, wn, wn, wn) x2 -> (B,) central voxel of each scaled target window."""
        em = torch.as_tensor(target_windows, device=self.device, dtype=self.dtype)
        mod = torch.as_tensor(reference_windows, device=self.device, dtype=self.dtype)
        B = em.shape[0]

        Fem = torch.fft.rfftn(em, dim=(1, 2, 3), norm="ortho").reshape(B, -1)
        Fmod = torch.fft.rfftn(mod, dim=(1, 2, 3), norm="ortho").reshape(B, -1)

        em_prof = self._radial_profile(Fem.abs())
        mod_prof = self._radial_profile(Fmod.abs())
        scale = mod_prof / em_prof
        # LocScale: scale_factor[~np.isfinite(scale_factor)] = 0 (empty/solvent windows)
        scale = torch.where(torch.isfinite(scale), scale, torch.zeros_like(scale))

        G = torch.zeros(B, self.n_prof, device=self.device, dtype=self.cdtype)
        G.index_add_(1, self._idx0, Fem * self._Wa)
        G.index_add_(1, self._idx1, Fem * self._Wb)
        return (scale.to(self.cdtype) * G).sum(dim=1).real


def _gather_windows(padded, centres, wn):
    """(n, wn, wn, wn) stack of windows whose centres sit at `centres` in the unpadded map."""
    out = np.empty((len(centres), wn, wn, wn), dtype=np.float32)
    for i, (k, j, l) in enumerate(centres):
        out[i] = padded[k:k + wn, j:j + wn, l:l + wn]
    return out


def run_windowed_scaling(target_map, reference_map, mask, wn=25, device="cpu",
                         chunk=4096, dtype=torch.float32, progress_callback=None):
    """Scale `target_map`'s local amplitudes to match `reference_map`, inside `mask`.

    Returns a map of the same shape: scaled values inside the mask, zero outside.

    `chunk` bounds memory -- roughly 344 KB per window at wn=25/float32, so 4096 windows is
    about 1.4 GB. Lower it if the GPU is small.
    """
    target_map = np.asarray(target_map, dtype=np.float32)
    reference_map = np.asarray(reference_map, dtype=np.float32)
    if target_map.shape != reference_map.shape:
        raise ValueError(f"map shapes differ: {target_map.shape} vs {reference_map.shape}")

    half = wn // 2
    # pad so a window can be taken around any masked voxel, including at the edges
    pad = [(half, wn - half)] * 3
    tgt_pad = np.pad(target_map, pad, mode="constant")
    ref_pad = np.pad(reference_map, pad, mode="constant")

    centres = np.argwhere(np.asarray(mask) > 0.5)
    if len(centres) == 0:
        raise ValueError("the mask selects no voxels -- nothing to scale")

    scaler = WindowedScaler(wn, device=device, dtype=dtype)
    output = np.zeros(target_map.shape, dtype=np.float32)

    for start in range(0, len(centres), chunk):
        sel = centres[start:start + chunk]
        tw = _gather_windows(tgt_pad, sel, wn)
        rw = _gather_windows(ref_pad, sel, wn)
        vals = scaler.scale_windows(tw, rw).cpu().numpy()
        output[sel[:, 0], sel[:, 1], sel[:, 2]] = vals
        if progress_callback is not None:
            progress_callback(min(start + chunk, len(centres)), len(centres))

    return output
