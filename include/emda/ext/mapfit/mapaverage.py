"""
Author: "Rangana Warshamanage, Garib N. Murshudov"
MRC Laboratory of Molecular Biology
    
This software is released under the
Mozilla Public License, version 2.0; see LICENSE.
"""

from __future__ import absolute_import, division, print_function, unicode_literals
import numpy as np
import fcodes_fast
from emda import core
from emda.ext.mapfit import map_class, emfit_class


def get_resol(bin_fsc, res_arr):
    bin_fsc = bin_fsc[bin_fsc > 0.1]
    if len(bin_fsc) > 0:
        dist = np.sqrt((bin_fsc - 0.143) ** 2)
    return res_arr[np.argmin(dist)]


def get_bfac(res_arr, signal):
    from scipy import stats

    high_res = 5.5  # Angstrom. This is a hard cutoff.
    res_ub = np.sqrt((res_arr - high_res) ** 2)
    ub_res = np.argmin(res_ub)
    res_arr_trunc, signal_trunc = trunc_array(res_arr, signal)
    if len(signal_trunc) < ub_res:
        s = 1 / res_arr[:ub_res]
        slope, intercept, _, _, _ = stats.linregress(
            (s * s) / 4.0, np.log(signal[:ub_res])
        )
    else:
        s = 1 / res_arr_trunc
        slope, intercept, _, _, _ = stats.linregress(
            (s * s) / 4.0, np.log(signal_trunc)
        )
    bfac = slope
    print("B factor = ", bfac, "Intercept = ", intercept)
    return bfac, intercept


def estimate_scale_and_relativebfac(fo_lst, uc, resol):
    from emda.ext.mapfit import curve_fit_3

    params = curve_fit_3.main(fo_lst, uc, resol)
    scale, bfac_relative = params
    return scale, bfac_relative


def set_array(arr, thresh=0.0):
    set2zero = False
    i = -1
    for ival in arr:
        i = i + 1
        if ival <= thresh and not set2zero:
            set2zero = True
        if set2zero:
            arr[i] = 0.0
    return arr


def trunc_array(res_arr, signal):
    val_lst = []
    res_lst = []
    j = -1
    fsc_thresh = 0.7  # this is an arbitrary choice.
    for ival in signal:
        if ival < 1.0:
            if ival > fsc_thresh:
                j = j + 1
                val_lst.append(ival)
                res_lst.append(res_arr[j])
            if ival < fsc_thresh:
                break
    return (
        np.array(res_lst, dtype="float", copy=False),
        np.array(val_lst, dtype="float", copy=False),
    )


def fsc_between_static_and_transfomed_map(
    staticmap, movingmap, bin_idx, rm, t, cell, nbin
):
    from emda.ext.mapfit import utils

    nx, ny, nz = staticmap.shape
    st, _, _, _ = fcodes_fast.get_st(nx, ny, nz, t)
    frt_full = utils.get_FRS(cell, rm, movingmap * st)[:, :, :, 0]
    f1f2_fsc = core.fsc.anytwomaps_fsc_covariance(staticmap, frt_full, bin_idx, nbin)[0]
    return f1f2_fsc


def get_ibin(bin_fsc):
    bin_fsc = bin_fsc[bin_fsc > 0.1]
    # dist = np.sqrt((bin_fsc - 0.143)**2)
    dist = np.sqrt((bin_fsc - 0.3) ** 2)
    ibin = np.argmin(dist) + 1
    return ibin


def main(
    maplist,
    ncycles,
    t_init,
    theta_init,
    smax,
    masklist,
    fobj,
    interp,
    fit=True,
    Bf_arr=None,
):
    from emda.ext.mapfit import newsignal, run_fit, utils

    if len(maplist) < 4:
        print("Map averaging requires at least two sets of maps " "with half sets")
        exit()
    elif len(maplist) >= 4 and len(maplist) % 2 == 0:
        print("Map average enabled")
        fobj.write("Map average enabled\n")
    try:
        emmap1 = map_class.EmmapAverage(hfmap_list=maplist, mask_list=masklist)
    except NameError:
        emmap1 = map_class.EmmapAverage(hfmap_list=maplist)
    emmap1.load_maps(fobj)
    emmap1.calc_fsc_variance_from_halfdata(fobj)
    # converting theta_init to rotmat for initial iteration
    fobj.write("\n")
    fobj.write("Initial fitting parameters:\n")
    fobj.write("    Translation: " + str(t_init) + " \n")
    fobj.write("    Rotation: " + str(theta_init) + " \n")
    q = core.quaternions.get_quaternion(theta_init)
    q = q / np.sqrt(np.dot(q, q))
    print("Initial quaternion: ", q)
    rotmat = core.quaternions.get_RM(q)
    fobj.write("    Rotation matrix: " + str(rotmat) + " \n")
    fobj.write("\n")
    fobj.write("    # fitting cycles: " + str(ncycles) + " \n")
    t = t_init
    if Bf_arr is None:
        Bf_arr=[0.0]
    rotmat_lst = []
    transl_lst = []
    # resolution estimate for line-fit
    dist = np.sqrt((emmap1.res_arr - smax) ** 2)
    slf = np.argmin(dist) + 1
    if slf % 2 != 0:
        slf = slf - 1
    slf = min([len(dist), slf])
    #
    if fit:
        for ifit in range(1, len(emmap1.eo_lst)):
            fobj.write("Fitting set#: " + str(ifit) + " \n")
            rotmat, t = run_fit.run_fit(
                emmap1, smax, rotmat, t, slf, ncycles, ifit, fobj, interp
            )
            rotmat_lst.append(rotmat)
            transl_lst.append(t)

        """newsignal.output_rotated_maps(emmap1,
                            rotmat_lst,
                            transl_lst,
                            fobj=fobj)"""
        """newsignal.output_rotated_maps_using_fo(emmap1,
                            rotmat_lst,
                            transl_lst,
                            fobj=fobj)"""
        newsignal.output_rotated_maps_using_fo_allmaps(
            emmap1, rotmat_lst, transl_lst, fobj=fobj
        )
    if not fit:
        """# using normalized Fourier coeffi
        cov_lst = []
        # estimating covaraince between current map vs. static map
        for frt in emmap1.fo_lst[1:]:
            bin_stats = \
                fsc.anytwomaps_fsc_covariance(emmap1.fo_lst[0], \
                    frt,emmap1.bin_idx, emmap1.nbin)
            f1f2_fsc, f1f2_covar = bin_stats[0], bin_stats[1]
            # mask f1f2_covar so that anything beyond zero get zeros
            f1f2_covar = set_array(f1f2_covar)
            cov_lst.append(f1f2_covar)        
        averagemaps = utils.avg_and_diffmaps(maps2avg=emmap1.eo_lst,
                                        uc=emmap1.map_unit_cell,
                                        nbin=emmap1.nbin,
                                        sgnl_var=emmap1.signalvar_lst,
                                        totl_var=emmap1.totalvar_lst,
                                        covar=cov_lst,
                                        hffsc=emmap1.hffsc_lst,
                                        bin_idx=emmap1.bin_idx,
                                        s_grid=emmap1.s_grid,
                                        res_arr=emmap1.res_arr,
                                        Bf_arr=Bf_arr)"""
        # using un-normalized Fourier coeffi.
        averagemaps = newsignal.avgmaps_using_fo_nofit(emmap1=emmap1, fobj=fobj)
        # averagemap output
        utils.output_maps(
            averagemaps=averagemaps,
            com_lst=[],
            t_lst=[],
            r_lst=[],
            unit_cell=emmap1.map_unit_cell,
            map_origin=emmap1.map_origin,
            bf_arr=Bf_arr,
            center=None,
        )
