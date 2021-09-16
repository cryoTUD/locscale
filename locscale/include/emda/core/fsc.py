"""
Author: "Rangana Warshamanage, Garib N. Murshudov"
MRC Laboratory of Molecular Biology
    
This software is released under the
Mozilla Public License, version 2.0; see LICENSE.
"""

from __future__ import absolute_import, division, print_function, unicode_literals

import numpy as np
import fcodes_fast
from emda.config import *


def halfmaps_fsc_variance(hf1, hf2, bin_idx, nbin):
    assert hf1.shape == hf2.shape
    nx, ny, nz = hf1.shape
    (
        fo,
        eo,
        noisevar,
        signalvar,
        totalvar,
        bin_fsc,
        bincount,
    ) = fcodes_fast.calc_fsc_using_halfmaps(
        hf1, hf2, bin_idx, nbin, debug_mode, nx, ny, nz
    )
    """ fo, eo, binstats, bincount = fcodes_fast.calc_halffsc(
        hf1, hf2, bin_idx, nbin, debug_mode, nx, ny, nz
    )
    noisevar = binstats[:,0]
    signalvar = binstats[:,1]
    totalvar = binstats[:,2]
    bin_fsc = binstats[:,3] """
    return bin_fsc, noisevar, signalvar, totalvar, fo, eo


def anytwomaps_fsc_covariance(f1, f2, bin_idx, nbin):
    assert f1.shape == f2.shape
    nx, ny, nz = f1.shape
    binstats, bin_arr_count = fcodes_fast.calc_fsc(
        f1, f2, bin_idx, nbin, debug_mode, nx, ny, nz
    )
    f1f2_covar = binstats[:, 0]
    bin_fsc = binstats[:, 1]

    """ f1f2_covar, bin_fsc, bincount = fcodes_fast.calc_covar_and_fsc_betwn_anytwomaps(
        f1, f2, bin_idx, nbin, debug_mode, nx, ny, nz
    ) """
    return [bin_fsc, f1f2_covar, bin_arr_count]


def predict_fsc(hf1, hf2, bin_idx, nbin, res_arr, nparticles=None, bfac=None):
    fsc, noise, signal, total, _, _ = halfmaps_fsc_variance(hf1, hf2, bin_idx, nbin)
    predicted_fsc_lst = []
    if nparticles is not None:
        factors = np.asarray(nparticles, dtype="float")
        print(factors)
        for i in range(len(factors)):
            predicted_fsc_lst.append(signal / (signal + noise * factors[i]))
        predicted_fsc_lst.append(2 * fsc / (1.0 + fsc))
        return predicted_fsc_lst
    if bfac is not None:
        bfac = np.asarray(bfac, dtype="float")
        print("B factors: ", bfac)
        s = ((1.0 / res_arr) ** 2) / 4.0
        for i in range(len(bfac)):
            sgnl_pred = signal * np.exp(bfac[i] * s)
            predicted_fsc_lst.append(sgnl_pred / (sgnl_pred + noise))
        predicted_fsc_lst.append(2 * fsc / (1.0 + fsc))
        return predicted_fsc_lst
