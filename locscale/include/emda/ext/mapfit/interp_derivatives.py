"""
Author: "Rangana Warshamanage, Garib N. Murshudov"
MRC Laboratory of Molecular Biology

This software is released under the
Mozilla Public License, version 2.0; see LICENSE.
"""

from __future__ import absolute_import, division, print_function, unicode_literals
import numpy as np
from emda import ext

# derivative calculation using interpolation
def interp_derivatives(dfs, rm):
    dFRs = ext.mapfit.utils.get_interp(rm, dfs[:, :, :, :], interp="cubic")
    return dFRs


def dfs_fullmap(mapin, xyz, vol):
    nx, ny, nz = mapin.shape
    dfs_full = np.zeros(shape=(nx, ny, nz, 3), dtype=np.complex64)
    for i in range(3):
        dfs_full[:, :, :, i] = np.fft.fftshift(
            (1 / vol) * 2j * np.pi * np.fft.fftn(mapin * xyz[i])
        )
    return dfs_full


def cut_dfs4interp(dfs_full, cx):
    nx = dfs_full.shape[0]
    dx = int((nx - 2 * cx) / 2)
    dy = int((nx - 2 * cx) / 2)
    dz = int((nx - 2 * cx) / 2)
    dfs = dfs_full[dx : dx + 2 * cx, dy : dy + 2 * cx, dz : dz + 2 * cx, :]
    return dfs
