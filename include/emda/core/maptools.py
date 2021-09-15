"""
Author: "Rangana Warshamanage, Garib N. Murshudov"
MRC Laboratory of Molecular Biology
    
This software is released under the
Mozilla Public License, version 2.0; see LICENSE.
"""

from __future__ import absolute_import, division, print_function, unicode_literals

import numpy as np
from emda.core import iotools, restools, plotter, fsc
import fcodes_fast
from emda.config import debug_mode


def test():
    print("maptools test ... Passed")


def Bf_linefit(res_arr, binSgnlVar, low_res_cutoff, high_res_cutoff, scale=1):
    from scipy import stats
    import numpy.ma as ma

    prediSigVar = np.zeros(len(res_arr), dtype="float")
    dist_low = np.sqrt((res_arr - low_res_cutoff) ** 2)
    lb = np.argmin(dist_low) + 1
    dist_high = np.sqrt((res_arr - high_res_cutoff) ** 2)
    ub = np.argmin(dist_high) + 1
    res_arr_ma = ma.masked_equal(res_arr, 0.0)
    s = 1 / res_arr_ma
    x_full = (s * s) / 2.0
    x_full = x_full.filled(0.0)
    slope, intercept, _, _, _ = stats.linregress(
        x_full[lb:ub], np.log(binSgnlVar[lb:ub])
    )
    Bfac = slope
    ln_s = -1.0 * abs(Bfac) * x_full + intercept
    prediSigVar = scale * np.exp(ln_s)
    return Bfac, prediSigVar


def estimate_map_resol(hfmap1, hfmap2):
    uc, arr1, _ = iotools.read_map(hfmap1)
    uc, arr2, _ = iotools.read_map(hfmap2)
    hf1 = np.fft.fftshift(np.fft.fftn(arr1))
    hf2 = np.fft.fftshift(np.fft.fftn(arr2))
    nbin, res_arr, bin_idx = restools.get_resolution_array(uc, hf1)
    bin_fsc, _, _, _, _, _ = fsc.halfmaps_fsc_variance(hf1, hf2, bin_idx, nbin)
    bin_fsc = bin_fsc[bin_fsc > 0.1]
    dist = np.sqrt((bin_fsc - 0.143) ** 2)
    map_resol = res_arr[np.argmin(dist)]
    return map_resol


def get_map_power(mapin):
    uc, arr, _ = iotools.read_map(mapin)
    hf = np.fft.fftshift(np.fft.fftn(arr))
    nx, ny, nz = hf.shape
    nbin, res_arr, bin_idx = restools.get_resolution_array(uc, hf)
    power_spectrum = fcodes_fast.calc_power_spectrum(
        hf, bin_idx, nbin, debug_mode, nx, ny, nz
    )
    print("Resolution   bin     Power")
    for i in range(len(res_arr)):
        print("{:.2f} {:.4f}".format(res_arr[i], power_spectrum[i]))
    return res_arr, power_spectrum


def get_biso_from_model(mmcif_file):
    _, _, _, _, Biso_np = iotools.read_mmcif(mmcif_file)
    Biso = np.median(Biso_np)
    return Biso


def get_biso_from_map(halfmap1, halfmap2):
    _, sig_var, res_arr = estimate_map_resol(halfmap1, halfmap2)
    plotter.plot_nlines_log(res_arr, [sig_var], ["Signal Variance"])
    low_res_cutoff = float(input("Enter low resolution cutoff: "))
    high_res_cutoff = float(input("Enter high resolution cutoff: "))
    Biso, preditcted_signal = Bf_linefit(
        res_arr, sig_var, low_res_cutoff, high_res_cutoff
    )
    plotter.plot_nlines_log(
        res_arr,
        [sig_var, preditcted_signal],
        ["Signal Variance", "Predicted SV"],
        "Predicted.eps",
    )
    return Biso


def apply_bfactor_to_map(mapname, bf_arr, mapout):
    uc, ar1, origin = iotools.read_map(mapname)
    hf1 = np.fft.fftshift(np.fft.fftn(ar1))
    nx, ny, nz = ar1.shape
    nbf = len(bf_arr)
    all_mapout = fcodes_fast.apply_bfactor_to_map(
        hf1, bf_arr, uc, debug_mode, nx, ny, nz, nbf
    )
    if mapout:
        for i in range(all_mapout.shape[3]):
            if bf_arr[i] < 0.0:
                Bcode = "_blur" + str(abs(bf_arr[i]))
            elif bf_arr[i] > 0.0:
                Bcode = "_sharp" + str(abs(bf_arr[i]))
            else:
                Bcode = "_unsharpened"
            filename_mrc = mapname[:-4] + Bcode + ".mrc"
            data2write = np.real(np.fft.ifftn(np.fft.ifftshift(all_mapout[:, :, :, i])))
            iotools.write_mrc(data2write, filename_mrc, uc, origin)
            print("writing done")
    return all_mapout


def map2mtz(mapname, mtzname="map2mtz.mtz", factor=1.0, resol=None):
    uc, ar1, _ = iotools.read_map(mapname)
    uc[0] = uc[0] * factor
    uc[1] = uc[1] * factor
    uc[2] = uc[2] * factor
    hf1 = np.fft.fftshift(np.fft.fftn(ar1))
    iotools.write_3d2mtz(unit_cell=uc, mapdata=hf1, outfile=mtzname, resol=resol)


def mtz2map(mtzname, map_size):
    from emda.core.mtz import mtz2map

    arr, unit_cell = mtz2map(mtzname=mtzname, map_size=map_size)
    return arr


def normalise_fo(fo, bin_idx=None, nbin=None, uc=None):
    from emda.core.restools import get_resolution_array
    from emda.core.fsc import halfmaps_fsc_variance

    if bin_idx is None:
        # calculate resolution grid
        nbin, res_arr, bin_idx = get_resolution_array(uc=uc, hf1=fo)
    # normalise Fo
    eo = halfmaps_fsc_variance(hf1=fo, hf2=fo, bin_idx=bin_idx, nbin=nbin)[5]
    return eo

