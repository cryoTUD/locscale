"""
Author: "Rangana Warshamanage, Garib N. Murshudov"
MRC Laboratory of Molecular Biology
    
This software is released under the
Mozilla Public License, version 2.0; see LICENSE.
"""

from __future__ import absolute_import, division, print_function, unicode_literals
from emda import core, ext
import numpy as np
import fcodes_fast
from emda.config import *


def lowpass_map(uc, arr, resol):
    # ideal filter
    fc = np.fft.fftshift(np.fft.fftn(arr))
    nbin, res_arr, bin_idx = core.restools.get_resolution_array(uc, fc)
    dist = np.sqrt((res_arr - resol) ** 2)
    cbin = np.argmin(dist) + 1
    fout = core.restools.cut_resolution(fc, bin_idx, res_arr, cbin)
    map_lwp = np.real(np.fft.ifftn(np.fft.ifftshift(fout)))
    return fout, map_lwp


def butterworth(uc, arr, smax, order=4):
    nx, ny, nz = arr.shape
    maxbin = np.amax(np.array([nx // 2, ny // 2, nz // 2]))
    nbin, res_arr, bin_idx, sgrid = fcodes_fast.resolution_grid(
        uc, debug_mode, maxbin, nx, ny, nz
    )
    fmap = np.fft.fftshift(np.fft.fftn(arr))
    dist = np.sqrt((res_arr[:nbin] - smax) ** 2)
    cbin = np.argmin(dist) + 1
    cutoffres = res_arr[cbin]
    order = 4  # order of the butterworth filter
    B = 1.0
    D = sgrid
    d = 1.0 / cutoffres
    bwfilter = 1.0 / (
        1 + B * ((D / d) ** (2 * order))
    )  # calculate the butterworth filter
    fmap_filtered = fmap * bwfilter
    map_lwp = np.real(np.fft.ifftn(np.fft.ifftshift(fmap_filtered)))
    return fmap_filtered, map_lwp
