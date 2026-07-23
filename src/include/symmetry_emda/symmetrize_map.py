# symmetrize map by operators
import numpy as np
import torch
from numpy.fft import fftn, ifftn, fftshift, ifftshift
from .GenerateOperators_v9_ky4 import operators_from_symbol, AngleAxis2rotatin
from .trilinear_torch import compute_nbin, pick_device, rotate_ft
"""
Original authors

Author: "Rangana Warshamanage, Garib N. Murshudov"
MRC Laboratory of Molecular Biology

https://gitlab.com/ccpem/emda/-/tree/master/
EMDA version 1.1.3.post6

The fcodes_fast Fortran kernels (resol_grid_em, trilinear2) are replaced by the
PyTorch implementation in trilinear_torch.py, validated against the compiled
originals to ~5e-15. Dropping the f2py/numpy.distutils build step lifts
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


def symmetrize_map(emmap, apix, pg="C1", twist=None, rise=None, n_steps=None, device=None, dtype=torch.complex128):
    """Symmetrise an in-memory map, dispatching on whether helical parameters are given.

    If `twist` is given, helical symmetry is imposed (combined with the Cn/Dn
    point group `pg`, see `symmetrize_map_helical`); `rise` is then required.
    Otherwise this is a plain point-group symmetrisation over `pg` (C1..Cn,
    D1..Dn, T, O, I), see `symmetrize_map_known_pg`.
    """
    if twist is not None:
        if rise is None:
            raise ValueError("`rise` must be provided together with `twist` for helical symmetrisation")
        return symmetrize_map_helical(emmap, apix, twist, rise, pg=pg, n_steps=n_steps, device=device, dtype=dtype)
    return symmetrize_map_known_pg(emmap, apix, pg, device=device, dtype=dtype)

def _parse_cyclic_or_dihedral(pg):
    """Validate that `pg` is a Cn or Dn symbol and return (order, is_dihedral).

    Helical symmetrisation only supports point groups whose principal axis can
    coincide with the helical axis (Cn) or whose dyad axes are perpendicular to
    it (Dn). T/O/I have no such axis and are not geometrically compatible
    with a helical rise/twist.
    """
    pg_l = pg.strip().upper()
    if len(pg_l) < 2 or pg_l[0] not in ("C", "D") or not pg_l[1:].isdigit():
        raise ValueError(
            "Helical symmetrisation only supports cyclic (Cn) or dihedral (Dn) "
            "point groups, got '{}'".format(pg)
        )
    order = int(pg_l[1:])
    if order < 1:
        raise ValueError("Point group order must be >= 1, got '{}'".format(pg))
    return order, pg_l[0] == "D"


def helical_operators(twist, rise, apix, n, pg="C1", n_steps=None):
    """Build (rotation, z-translation) pairs for combined helical + Cn/Dn symmetry.

    The helical axis is the volume's first (z) array axis, which is the same axis
    used by the xyz-convention point-group operators from
    ``operators_from_symbol`` once passed through ``_to_zyx``. A Cn point
    group is therefore coincident with the helical axis, and the dyad of a Dn
    point group is perpendicular to it, matching the standard convention for
    n-start helices (e.g. RELION/cryoSPARC helical symmetry model).

    Parameters
    ----------
    twist : float
        Helical twist in degrees per asymmetric unit. Sign gives handedness
        (following the cryoSPARC rise/twist convention); flip it if the
        result comes out the wrong hand.
    rise : float
        Helical rise in Angstrom per asymmetric unit (> 0).
    apix : float
        Isotropic voxel size in Angstrom/pixel.
    n : int
        Box size (voxels) along the axis the translation is applied on, i.e.
        the size of the (padded) volume the operators will be used on.
    pg : str
        "C<n>" or "D<n>" internal point-group symmetry of the helical
        assembly, coincident with (Cn) or perpendicular to (Dn dyad) the
        helical axis. Default "C1" (no extra symmetry beyond the helix).
    n_steps : int, optional
        Number of times the (twist, rise) screw operation is applied on each
        side of the reference copy (k = -n_steps .. +n_steps). This is a step
        count, not an asymmetric-unit count: with an internal Cn/Dn point
        group (n-start helix), each step still covers `pg`'s full order in
        asymmetric units, since every step is combined with every point-group
        operator. It's also not a "repeat" in the crystallographic sense (the
        axial distance at which twist exactly recurs mod 360 degrees) or the
        pitch (rise per 360-degree turn) -- those are generally different,
        sometimes non-finite, quantities. Defaults to the largest value that
        keeps every translated copy inside a zero-padded box of size n (i.e.
        no FFT wraparound); see `double_the_axes`/`rebox_map`.

    Returns
    -------
    list of (np.ndarray(3,3), float)
        Rotation (xyz convention, to be passed through `apply_op`/`_to_zyx`)
        and z-translation in pixels, for each combined symmetry operator.
    """
    if rise <= 0:
        raise ValueError("rise must be > 0, got {}".format(rise))
    _parse_cyclic_or_dihedral(pg)  # validates pg is Cn/Dn; raises otherwise
    _, _, pg_ops = operators_from_symbol(pg)

    rise_px = rise / apix
    if n_steps is None:
        # Half of `n` is the zero-padding introduced by double_the_axes; stay
        # inside it so translated copies never wrap around the FFT torus.
        n_steps = max(0, int(np.floor((n / 2.0) / rise_px)))

    twist_rad = np.deg2rad(twist)
    axis_z = np.array([0.0, 0.0, 1.0])

    ops = []
    for k in range(-n_steps, n_steps + 1):
        R_k = AngleAxis2rotatin(axis_z, k * twist_rad)
        t_k = np.array([0.0, 0.0, k * rise_px])
        for R_p in pg_ops:
            combined_R = R_p @ R_k
            combined_t = R_p @ t_k
            ops.append((combined_R, float(combined_t[2])))
    return ops


def apply_helical_op(f1, op, dz_pixels, nbin, device=None, dtype=torch.complex128):
    """Rotate (as `apply_op`) then apply the z-translation phase ramp (shift theorem)."""
    rotated = apply_op(f1, op, nbin, device=device, dtype=dtype)
    n = rotated.shape[0]
    dev = rotated.device
    if dz_pixels == 0:
        return rotated
    real_dtype = torch.float32 if dtype == torch.complex64 else torch.float64
    nmin = -(n // 2)
    idx = torch.arange(nmin, nmin + n, device=dev, dtype=real_dtype)
    angle = (-2.0 * np.pi * dz_pixels / n) * idx
    phase = torch.complex(torch.cos(angle), torch.sin(angle))
    return rotated * phase.view(-1, 1, 1)


def symmetrize_map_helical(emmap, apix, twist, rise, pg="C1", n_steps=None, device=None, dtype=torch.complex128):
    """Average a map over helical symmetry, optionally combined with a Cn/Dn point group.

    The helical axis is assumed to be the volume's first (z) array axis. The
    box is zero-padded (`double_the_axes`) before the z-translations are
    applied, to keep them from wrapping around the FFT torus, and cropped
    back (`rebox_map`) to the original size afterwards. Because of this finite
    averaging window, voxels near the top/bottom of the box receive fewer
    effective symmetry-related contributions than the centre. Provide a box
    taller than the region of interest if this matters.

    Parameters
    ----------
    emmap : (n, n, n) ndarray
    apix : float, Angstrom/pixel (isotropic)
    twist : float, helical twist in degrees per asymmetric unit
    rise : float, helical rise in Angstrom per asymmetric unit (> 0)
    pg : str, "C<n>" or "D<n>" internal point-group symmetry (default "C1")
    n_steps : int, optional, see `helical_operators`
    """
    print("===== Symmetrize Map (helical) =====")
    print("Credits: Rangana Warshamanage, Garib N. Murshudov")
    print("EMDA version 1.1.3.post6")
    print("https://gitlab.com/ccpem/emda/-/tree/master/")
    print("======================================")

    padded = double_the_axes(emmap)
    n = padded.shape[0]
    dev = pick_device(device)

    ops = helical_operators(twist, rise, apix, n, pg=pg, n_steps=n_steps)
    print("Symmetrising {} map (padded to {}) over {} helical x point-group operators on {}".format(
        emmap.shape, padded.shape, len(ops), dev))

    f1 = fftshift(fftn(fftshift(padded)))
    nbin = compute_nbin(n)

    f1_t = torch.as_tensor(f1, dtype=dtype, device=dev)
    frs_sum = torch.zeros_like(f1_t)
    for op, dz_pixels in ops:
        frs_sum += apply_helical_op(f1_t, op, dz_pixels, nbin, device=dev, dtype=dtype)
    avg_f = (frs_sum / len(ops)).cpu().numpy()

    avgmap = ifftshift(np.real(ifftn(ifftshift(avg_f))))
    avgmap = rebox_map(avgmap)
    return avgmap
