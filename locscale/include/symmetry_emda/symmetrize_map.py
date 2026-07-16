# symmetrize map by operators
import numpy as np
import torch
from numpy.fft import fftn, ifftn, fftshift, ifftshift
from locscale.include.symmetry_emda.GenerateOperators_v9_ky4 import operators_from_symbol
from locscale.include.symmetry_emda.trilinear_torch import compute_nbin, pick_device, rotate_ft
"""
Original authors

Author: "Rangana Warshamanage, Garib N. Murshudov"
MRC Laboratory of Molecular Biology

https://gitlab.com/ccpem/emda/-/tree/master/
EMDA version 1.1.3.post6

The fcodes_fast Fortran kernels (resol_grid_em, trilinear2) are replaced by the
PyTorch implementation in trilinear_torch.py, validated against the compiled
originals to ~5e-15. Dropping the f2py/numpy.distutils build step is what lifts
the numpy<=1.26 pin and allows a compiler-free install.
"""

def double_the_axes(arr1):
    nx, ny, nz = arr1.shape
    big_arr1 = np.zeros((2 * nx, 2 * ny, 2 * nz), dtype="float")
    dx = int(nx / 2)
    dy = int(ny / 2)
    dz = int(nz / 2)
    big_arr1[dx : dx + nx, dy : dy + ny, dz : dz + nz] = arr1
    return big_arr1


def _to_zyx(op):
    """Reverse an xyz-convention operator into the volume's zyx axis order.

    Reproduces the row/column reversal applied before the Fortran trilinear2
    call, i.e. rm = J @ op @ J with J the exchange matrix.
    """
    assert op.ndim == 2
    assert op.shape[0] == op.shape[1] == 3
    tmp = np.zeros(op.shape, 'float')
    rm = np.zeros(op.shape, 'float')
    tmp[:,0] = op[:,2]
    tmp[:,1] = op[:,1]
    tmp[:,2] = op[:,0]
    rm[0, :] = tmp[2, :]
    rm[1, :] = tmp[1, :]
    rm[2, :] = tmp[0, :]
    return rm


def apply_op(f1, op, nbin, device=None, dtype=torch.complex128):
    rm = _to_zyx(op)
    return rotate_ft(f1, rm, nbin=nbin, device=device, dtype=dtype)


def rebox_map(arr1):
    nx, ny, nz = arr1.shape
    dx = int(nx / 4)
    dy = int(ny / 4)
    dz = int(nz / 4)
    reboxed_map = arr1[dx : dx + nx//2, dy : dy + ny//2, dz : dz + nz//2]
    return reboxed_map


def symmetrize_map_known_pg(emmap, apix, pg, device=None, dtype=torch.complex128):
    print("===== Symmetrize Map =====")
    print("Credits: Rangana Warshamanage, Garib N. Murshudov")
    print("EMDA version 1.1.3.post6")
    print("https://gitlab.com/ccpem/emda/-/tree/master/")
    print("==========================")

    _, _, ops = operators_from_symbol(pg)
    dev = pick_device(device)
    print("Symmetrising {} map over {} operators on {}".format(emmap.shape, len(ops), dev))

    f1 = fftshift(fftn(fftshift(emmap)))
    nbin = compute_nbin(f1.shape[0])

    # Accumulate on-device; only the averaged volume returns to the host.
    f1_t = torch.as_tensor(f1, dtype=dtype, device=dev)
    frs_sum = torch.zeros_like(f1_t)
    for op in ops:
        frs_sum += apply_op(f1_t, op, nbin, device=dev, dtype=dtype)
    avg_f = (frs_sum / len(ops)).cpu().numpy()

    avgmap = ifftshift(np.real(ifftn(ifftshift(avg_f))))
    #avgmap = rebox_map(avgmap)
    return avgmap


def symmetrize_map_emda(emmap_path, pg, device=None, dtype=torch.complex128):
    from locscale.include.emmer.ndimage.map_utils import load_map
    emmap,apix = load_map(emmap_path)
    symmetry_average_map = symmetrize_map_known_pg(emmap, apix, pg, device=device, dtype=dtype)

    return symmetry_average_map
