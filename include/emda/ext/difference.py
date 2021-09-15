"""
Author: "Rangana Warshamanage, Garib N. Murshudov"
MRC Laboratory of Molecular Biology
    
This software is released under the
Mozilla Public License, version 2.0; see LICENSE.
"""

# difference map
from __future__ import absolute_import, division, print_function, unicode_literals
import numpy as np
import fcodes_fast
from emda import core, ext



def diffmap_scalebypower(
    f1,
    f2,
    cell,
    smax=0.0,
    origin=None,
    res_arr=None,
    bin_idx=None,
    nbin=None,
    s_grid=None,
):

    if bin_idx is None:
        nbin, res_arr, bin_idx = core.restools.get_resolution_array(cell, f1)

    scale = ext.scale_maps.scale_twomaps_by_power(
        f1=f1, f2=f2, bin_idx=bin_idx, res_arr=res_arr
    )
    if smax >= 0.0:
        mask = res_arr > smax
        for i in range(len(res_arr)):
            res = res_arr[i] * mask[i]
            scl = scale[i] * mask[i]
            print(res, scl)
    else:
        print("resolution cannot be negative. Stopping now...")
        exit()

    scale_grid = ext.scale_maps.transfer_power(bin_idx, res_arr, scale * mask)
    f2_scaled = scale_grid * f2
    mask_grid = ext.scale_maps.transfer_power(bin_idx, res_arr, mask)
    f1 = f1 * mask_grid
    # calculating maps
    dm1_dm2 = np.real(np.fft.ifftshift(np.fft.ifftn(np.fft.ifftshift(f1 - f2_scaled))))
    dm2_dm1 = np.real(np.fft.ifftshift(np.fft.ifftn(np.fft.ifftshift(f2_scaled - f1))))
    core.iotools.write_mrc(dm1_dm2, "diffmap_m1-m2_pwr.mrc", cell, origin)
    core.iotools.write_mrc(dm2_dm1, "diffmap_m2-m1_pwr.mrc", cell, origin)
    return dm1_dm2, dm2_dm1


def diffmap_scalebyampli(
    f1, f2, cell, smax=0.0, origin=None, res_arr=None, bin_idx=None, nbin=None
):

    if bin_idx is None:
        nbin, res_arr, bin_idx = core.restools.get_resolution_array(cell, f1)
    nx, ny, nz = f1.shape
    diffmap = fcodes_fast.differencemap(
        fo=f1,
        fc=f2,
        bin_idx=bin_idx,
        res_arr=res_arr,
        smax=smax,
        mode=1,
        nbin=len(res_arr),
        nx=nx,
        ny=ny,
        nz=nz,
    )
    return diffmap


def diffmap_normalisedsf(
    f1, f2, cell, smax=0.0, origin=None, res_arr=None, bin_idx=None, nbin=None
):

    if bin_idx is None:
        nbin, res_arr, bin_idx = core.restools.get_resolution_array(cell, f1)

    nx, ny, nz = f1.shape

    diffmap = fcodes_fast.diffmap_norm(
        fo=f1,
        fc=f2,
        bin_idx=bin_idx,
        res_arr=res_arr,
        smax=smax,
        mode=1,
        nbin=len(res_arr),
        nx=nx,
        ny=ny,
        nz=nz,
    )
    return diffmap
