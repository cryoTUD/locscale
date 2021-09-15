"""
Author: "Rangana Warshamanage, Garib N. Murshudov"
MRC Laboratory of Molecular Biology
    
This software is released under the
Mozilla Public License, version 2.0; see LICENSE.
"""

# curve fit
#Insights: 
#https://stackoverflow.com/questions/16745588/least-squares-minimization-complex-numbers/20104454#20104454
#Method of minimization: scipy.optimize.leastsq; Accoding to scipy documentation method os Levenberg–Marquardt

from __future__ import absolute_import, division, print_function, unicode_literals

import numpy as np
import fcodes_fast
from emda import core
from emda.core import restools, fsc, iotools


def model(x, f2, s):
    k = x[0]
    b = x[1]
    return k * f2 * np.exp(-b * s)


def fun(x, f1, f2, s):
    diff = model(x, f2, s) - f1
    z1d = np.zeros(f1.size * 2, dtype=np.float64)
    z1d[0 : z1d.size : 2] = diff.real
    z1d[1 : z1d.size : 2] = diff.imag
    return z1d


def jac(x, f1, f2, s):
    J = np.zeros((f1.size * 2, 2), dtype=np.float64)
    bs = x[1] * s
    ks = x[0] * s
    j0 = np.exp(-bs) * f2
    J[0 : f1.size * 2 : 2, 0] = j0.real
    J[1 : f1.size * 2 : 2, 0] = j0.imag
    j1 = -ks * np.exp(-bs) * f2
    J[0 : f1.size * 2 : 2, 1] = j1.real
    J[1 : f1.size * 2 : 2, 1] = j1.imag
    return J


def lsq(f1, f2, s, x0):
    from scipy.optimize import leastsq

    params, _, _, _, _ = leastsq(
        fun, x0, Dfun=jac, args=(f1, f2, s), full_output=True
    )
    return params


def get_resolution(fhf_lst, uc):
    assert fhf_lst[0].shape == fhf_lst[1].shape
    nbin, res_arr, bin_idx = restools.get_resolution_array(uc, fhf_lst[0])
    bin_fsc = fsc.anytwomaps_fsc_covariance(
        fhf_lst[0], fhf_lst[1], bin_idx, nbin
    )[0]
    bin_fsc = bin_fsc[bin_fsc > 0.1]
    if len(bin_fsc) > 0:
        dist = np.sqrt((bin_fsc - 0.143) ** 2)
        resol = res_arr[np.argmin(dist)]
        return resol
    else:
        raise SystemExit("Length of resolution array is zero!!!")


def main(fhf_lst, uc, resol=5.0):
    from scipy.optimize import leastsq

    assert fhf_lst[0].shape == fhf_lst[1].shape
    nx, ny, nz = fhf_lst[0].shape
    print(nx, ny, nz)
    maxbin = np.amax(np.array([nx // 2, ny // 2, nz // 2]))
    _, s_grid, mask = fcodes_fast.resolution_grid_full(uc, resol, 1, maxbin, nx, ny, nz)
    # optimization
    x0 = np.array([1.0, 10.0])
    f1 = (fhf_lst[0] * mask).flatten()
    f2 = (fhf_lst[1] * mask).flatten()
    s = s_grid.flatten()
    s = (s ** 2) / 4
    params, _, _, _, _ = leastsq(
        fun, x0, Dfun=jac, args=(f1, f2, s), full_output=True
    )
    return params


if __name__ == "__main__":
    maplist = [
        "/Users/ranganaw/MRC/Map_superposition/testmaps/str11.map",
        "/Users/ranganaw/MRC/Map_superposition/testmaps/emda/hf1_rotated_5deg.mrc",
    ]
    #'/Users/ranganaw/MRC/Map_superposition/testmaps/apply_bfac/avgmap__sharp50.mrc']
    fhf_lst = []
    for imap in maplist:
        uc, arr, _ = iotools.read_map(imap)
        fhf_lst.append(np.fft.fftshift(np.fft.fftn(np.fft.fftshift(arr))))
    resol = get_resolution(fhf_lst, uc)
    main(fhf_lst, uc, resol)

