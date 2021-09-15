"""
Author: "Rangana Warshamanage, Garib N. Murshudov"
MRC Laboratory of Molecular Biology
    
This software is released under the
Mozilla Public License, version 2.0; see LICENSE.
"""

from __future__ import absolute_import, division, print_function, unicode_literals
from emda.core import restools


def frequency_marching(fo, bin_idx, res_arr, bmax=None, fobj=None):
    cbin = cx = bmax
    print("cbin=", cbin)
    print("fit resolution:", res_arr[cbin])
    if fobj is not None:
        fobj.write("fit resolution: " + str(res_arr[cbin]) + " (A) \n")
    print()
    dx = int((fo.shape[0] - 2 * cx) / 2)
    dy = int((fo.shape[1] - 2 * cx) / 2)
    dz = int((fo.shape[2] - 2 * cx) / 2)
    cBIdx = bin_idx[dx : dx + 2 * cx, dy : dy + 2 * cx, dz : dz + 2 * cx]
    cutmap = restools.cut_resolution(fo, bin_idx, res_arr, cbin)[
        dx : dx + 2 * cx, dy : dy + 2 * cx, dz : dz + 2 * cx
    ]
    return cutmap, cBIdx, cbin

