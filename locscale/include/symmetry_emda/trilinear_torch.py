"""PyTorch replacement for the fcodes_fast (EMDA) Fortran routines used in symmetrisation.

This replaces two f2py entry points:

  * ``resol_grid_em``  -> :func:`compute_nbin` + the radius test in :func:`rotate_ft`
  * ``trilinear2``     -> :func:`rotate_ft`

Removing the Fortran is what lets LocScale build as a pure-Python wheel (no
compiler, no numpy.distutils) and therefore install on numpy>=2 and on Colab.

Fidelity: validated against the compiled fcodes_fast on cubic volumes at
n = 16/32/48 over random SO(3) rotations; worst relative deviation ~5e-15,
i.e. float64 round-off. See tests/test_symmetry_torch.py.

Conventions inherited from the Fortran, which the caller depends on:
  * The volume is centered (zero frequency at index n//2), as produced by
    ``fftshift(fftn(fftshift(map)))``; index range is [-n//2, (n-2)//2].
  * Sampling is a pull-back through ``rm^T`` (Fortran: matmul(transpose(RM), s)).
  * Only the half-space h <= 0 is interpolated; the rest is filled by Friedel
    symmetry, skipping the Nyquist planes that have no partner inside the box.
"""

import numpy as np
import torch


def compute_nbin(n):
    """Number of resolution shells for a cubic box of size `n`.

    Mirrors resol_grid_em: xyzmax = n//2 - 1 and the shell loop runs
    i = 0 .. xyzmax-1, so nbin = xyzmax.
    """
    return n // 2 - 1


def pick_device(device=None):
    """Resolve the compute device, preferring CUDA when available.

    MPS is never auto-selected: Apple GPUs have no float64, so the default
    complex128 kernel cannot run there. Pass device="mps" with
    dtype=torch.complex64 to opt in.
    """
    if device is not None:
        return torch.device(device)
    if torch.cuda.is_available():
        return torch.device("cuda")
    return torch.device("cpu")


def rotate_ft(F, rm, nbin=None, device=None, chunk=None, dtype=torch.complex128):
    """Rotate a centered Hermitian Fourier volume by `rm` with trilinear interpolation.

    Parameters
    ----------
    F : (n, n, n) complex ndarray or tensor
        Centered Fourier volume.
    rm : (3, 3) array
        Rotation matrix in the same axis convention the Fortran received, i.e.
        already reversed to (z, y, x) by the caller.
    nbin : int, optional
        Shell count; defaults to :func:`compute_nbin`.
    device : str or torch.device, optional
        Defaults to CUDA when available, else CPU.
    chunk : int, optional
        Number of h-planes processed at once. Bounds peak memory, which would
        otherwise be 8 gathered complex volumes at once (~1 GB at n=256).
    dtype : torch.dtype, optional
        complex128 (default) reproduces the Fortran to round-off. complex64 is
        markedly faster on GPUs, which run float64 at a fraction of float32 rate,
        and is required on MPS; it costs ~1e-7 relative accuracy.

    Returns
    -------
    torch.Tensor : (n, n, n) complex tensor of `dtype` on `device`.
    """
    dev = pick_device(device)
    if dev.type == "mps" and dtype == torch.complex128:
        raise ValueError(
            "MPS does not support float64/complex128; pass dtype=torch.complex64 "
            "to run on Apple GPUs, or use device='cpu' for full precision."
        )
    real_dtype = torch.float32 if dtype == torch.complex64 else torch.float64
    Ft = torch.as_tensor(np.asarray(F) if not torch.is_tensor(F) else F,
                         dtype=dtype, device=dev)
    if Ft.ndim != 3 or not (Ft.shape[0] == Ft.shape[1] == Ft.shape[2]):
        raise ValueError(f"expected a cubic 3D volume, got shape {tuple(Ft.shape)}")

    n = Ft.shape[0]
    nmin, nmax = -(n // 2), (n - 2) // 2
    if nbin is None:
        nbin = compute_nbin(n)
    rmt = torch.as_tensor(np.asarray(rm), dtype=real_dtype, device=dev)

    FRS = torch.zeros_like(Ft)
    radius_max = (nbin - 1) + 2.5

    h_all = torch.arange(nmin, 1, device=dev, dtype=real_dtype)
    line = torch.arange(nmin, nmax + 1, device=dev, dtype=real_dtype)
    if chunk is None:
        # ~64 MiB per gathered corner at complex128.
        chunk = max(1, int(4_000_000 / (n * n)))

    for start in range(0, h_all.numel(), chunk):
        h = h_all[start:start + chunk]
        H, K, L = torch.meshgrid(h, line, line, indexing="ij")
        S = torch.stack([H, K, L], dim=-1)

        valid = torch.sqrt((S ** 2).sum(-1)) <= radius_max

        # x_i = sum_j rm[j, i] * s_j   (Fortran: matmul(transpose(RM), s))
        X = torch.einsum("ji,...j->...i", rmt, S)
        X0 = torch.floor(X)
        X1 = X0 + 1.0
        # Both corners must be inside the box; otherwise the point stays zero.
        valid &= ((X0 >= nmin) & (X0 <= nmax) & (X1 >= nmin) & (X1 <= nmax)).all(-1)

        XD = X - X0
        A = (X0 - nmin).long().clamp_(0, n - 1)
        B = (X1 - nmin).long().clamp_(0, n - 1)
        a0, a1, a2 = A[..., 0], A[..., 1], A[..., 2]
        b0, b1, b2 = B[..., 0], B[..., 1], B[..., 2]

        w0 = XD[..., 0].to(Ft.dtype); v0 = (1.0 - XD[..., 0]).to(Ft.dtype)
        w1 = XD[..., 1].to(Ft.dtype); v1 = (1.0 - XD[..., 1]).to(Ft.dtype)
        w2 = XD[..., 2].to(Ft.dtype); v2 = (1.0 - XD[..., 2]).to(Ft.dtype)

        c00 = Ft[a0, a1, a2] * v0 + Ft[b0, a1, a2] * w0
        c01 = Ft[a0, a1, b2] * v0 + Ft[b0, a1, b2] * w0
        c10 = Ft[a0, b1, a2] * v0 + Ft[b0, b1, a2] * w0
        c11 = Ft[a0, b1, b2] * v0 + Ft[b0, b1, b2] * w0
        c = (c00 * v1 + c10 * w1) * v2 + (c01 * v1 + c11 * w1) * w2
        c = torch.where(valid, c, torch.zeros((), dtype=Ft.dtype, device=dev))

        FRS[(H - nmin).long(), (K - nmin).long(), (L - nmin).long()] = c

        m = valid & (H != nmin) & (K != nmin) & (L != nmin)
        FRS[(-H[m] - nmin).long(), (-K[m] - nmin).long(), (-L[m] - nmin).long()] = torch.conj(c[m])

    return FRS
