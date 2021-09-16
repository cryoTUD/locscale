"""
Author: "Rangana Warshamanage, Garib N. Murshudov"
MRC Laboratory of Molecular Biology

This software is released under the
Mozilla Public License, version 2.0; see LICENSE.
"""

# Applying B-factor and calculating signal

from __future__ import absolute_import, division, print_function, unicode_literals
import numpy as np
import fcodes_fast
from emda import core
from emda.ext.mapfit import mapaverage, utils, run_fit


def output_rotated_maps(emmap1, r_lst, t_lst, Bf_arr=None, fobj=None):
    from scipy.ndimage.interpolation import shift

    if Bf_arr is None:
        Bf_arr = [0.0]
    com_lst = emmap1.com_lst
    fo_lst = emmap1.fo_lst
    eo_lst = emmap1.eo_lst
    hffsc_lst = []
    hffsc_lst.append(mapaverage.set_array(emmap1.hffsc_lst[0]))
    totalvar_lst = []
    totalvar_lst.append(emmap1.totalvar_lst[0])
    bin_idx = emmap1.bin_idx
    nbin = emmap1.nbin
    res_arr = emmap1.res_arr
    cell = emmap1.map_unit_cell
    origin = emmap1.map_origin
    nx, ny, nz = fo_lst[0].shape
    frt_lst = []
    frt_lst.append(fo_lst[0])
    ert_lst = []
    ert_lst.append(eo_lst[0])
    cov_lst = []
    signalvar_lst = []
    # static map
    data2write = np.real(np.fft.ifftshift(np.fft.ifftn(np.fft.ifftshift(fo_lst[0]))))
    data2write = shift(data2write, np.subtract(com_lst[0], emmap1.box_centr))
    core.iotools.write_mrc(data2write, "static_map.mrc", cell, origin)
    fobj.write("\n static_map.mrc was written. \n")
    signalvar_static = emmap1.signalvar_lst[0]
    fobj.write("\n ***Map averaging using normalized structure factors*** \n")
    # smoothening signal
    fobj.write("\nStatic map signal \n")
    signalvar_static = get_extended_signal(
        res_arr, mapaverage.set_array(signalvar_static), hffsc_lst[0], fobj=fobj
    )
    signalvar_lst.append(mapaverage.set_array(signalvar_static))
    i = -1
    for t, rotmat in zip(t_lst, r_lst):
        i = i + 1
        st, _, _, _ = fcodes_fast.get_st(nx, ny, nz, t)
        # apply transformation on half maps of moving map
        fobj.write("Apply transformation on half maps... \n")
        frt_hf1 = utils.get_FRS(cell, rotmat, emmap1.fhf_lst[2 * i + 2] * st)[
            :, :, :, 0
        ]
        frt_hf2 = utils.get_FRS(cell, rotmat, emmap1.fhf_lst[2 * i + 3] * st)[
            :, :, :, 0
        ]
        frt_full = (frt_hf1 + frt_hf2) / 2.0
        data2write = np.real(np.fft.ifftshift(np.fft.ifftn(np.fft.ifftshift(frt_full))))
        data2write = shift(data2write, np.subtract(com_lst[0], emmap1.box_centr))
        core.iotools.write_mrc(
            data2write, "{0}_{1}.{2}".format("fitted_map", str(i), "mrc"), cell, origin
        )
        # estimating covariance between current map vs. static map
        _, f1f2_fsc = fcodes_fast.calc_covar_and_fsc_betwn_anytwomaps(
            fo_lst[0], frt_full, bin_idx, nbin, 0, nx, ny, nz
        )
        # find the relative B-factor between static and frt_full
        resol = mapaverage.get_resol(f1f2_fsc, res_arr)
        scale, bfac_relative = mapaverage.estimate_scale_and_relativebfac(
            [frt_lst[0], frt_full], cell, resol
        )
        fobj.write(
            "Relative B-factor between static map and moving(full) map: "
            + str(bfac_relative)
            + "\n"
        )
        # apply B-factor and scale to 2nd map
        s_grid = (emmap1.s_grid ** 2) / 4.0
        frt_full_scaled = scale * frt_full * np.exp(bfac_relative * s_grid)
        # re-estimate covariance
        f1f2_covar, f1f2_fsc = fcodes_fast.calc_covar_and_fsc_betwn_anytwomaps(
            fo_lst[0], frt_full_scaled, bin_idx, nbin, 0, nx, ny, nz
        )
        # smoothen signal
        fobj.write("\nCovariance between static and moving maps \n")
        f1f2_covar = get_extended_signal(
            res_arr, mapaverage.set_array(f1f2_covar), f1f2_fsc, fobj=fobj
        )
        # mask f1f2_covar so that anything beyond zero get zeros
        cov_lst.append(mapaverage.set_array(f1f2_covar))
        # fsc12_lst.append(f1f2_fsc)
        scale1, bfac_relative1 = mapaverage.estimate_scale_and_relativebfac(
            [frt_lst[0], frt_hf1], cell, resol
        )
        scale2, bfac_relative2 = mapaverage.estimate_scale_and_relativebfac(
            [frt_lst[0], frt_hf2], cell, resol
        )
        frt_hf1_scaled = scale1 * frt_hf1 * np.exp(bfac_relative1 * s_grid)
        frt_hf2_scaled = scale2 * frt_hf2 * np.exp(bfac_relative2 * s_grid)
        (
            bin_fsc,
            _,
            signalvar,
            totalvar,
            frt_full_scaled,
            ert_full,
        ) = core.fsc.halfmaps_fsc_variance(
            frt_hf1_scaled, frt_hf2_scaled, bin_idx, nbin
        )
        ert_lst.append(ert_full)
        # smoothen signal
        bin_fsc_full = 2 * bin_fsc / (1 + bin_fsc)
        fobj.write("\nMoving map signal \n")
        signalvar = get_extended_signal(
            res_arr, mapaverage.set_array(signalvar), bin_fsc_full, fobj=fobj
        )
        signalvar_lst.append(mapaverage.set_array(signalvar))
        totalvar_lst.append(mapaverage.set_array(totalvar))
        hffsc_lst.append(mapaverage.set_array(bin_fsc_full))

    var_cov_lst = signalvar_lst + cov_lst
    core.plotter.plot_nlines_log(res_arr, var_cov_lst)
    fvarcov = open("var_cov.txt", "w")
    for i in range(nbin):
        s11 = var_cov_lst[0][i]
        s22 = var_cov_lst[1][i]
        s12 = var_cov_lst[2][i]
        fvarcov.write("{:.2f} {:.4f} {:.4f} {:.4f}\n".format(res_arr[i], s11, s22, s12))
    # map averaging with normalised maps
    averagemaps = utils.avg_and_diffmaps(
        ert_lst,
        cell,
        nbin,
        signalvar_lst,
        totalvar_lst,
        cov_lst,
        hffsc_lst,
        bin_idx,
        emmap1.s_grid,
        res_arr,
        Bf_arr,
    )

    utils.output_maps(
        averagemaps=averagemaps,
        com_lst=emmap1.com_lst,
        t_lst=t_lst,
        r_lst=r_lst,
        center=emmap1.box_centr,
        unit_cell=cell,
        map_origin=origin,
        bf_arr=Bf_arr,
    )


def output_rotated_maps_using_fo(emmap1, r_lst, t_lst, Bf_arr=None, fobj=None):
    from scipy.ndimage.interpolation import shift

    if Bf_arr is None:
        Bf_arr = [0.0]
    com_lst = emmap1.com_lst
    fo_lst = emmap1.fo_lst
    signalvar_lst = []
    totalvar_lst = []
    totalvar_lst.append(emmap1.totalvar_lst[0])
    hffsc_lst = []
    hffsc_lst.append(mapaverage.set_array(emmap1.hffsc_lst[0]))
    bin_idx = emmap1.bin_idx
    nbin = emmap1.nbin
    res_arr = emmap1.res_arr
    cell = emmap1.map_unit_cell
    origin = emmap1.map_origin
    nx, ny, nz = fo_lst[0].shape
    frt_lst = []
    frt_lst.append(fo_lst[0])
    cov_lst = []
    # static map
    data2write = np.real(np.fft.ifftshift(np.fft.ifftn(np.fft.ifftshift(fo_lst[0]))))
    data2write = shift(data2write, np.subtract(com_lst[0], emmap1.box_centr))
    core.iotools.write_mrc(data2write, "static_map.mrc", cell, origin)
    del data2write
    fobj.write("\n *** Map averaging using un-normalized structure factors *** \n")
    maps2average_lst = []
    maps2average_lst.append(fo_lst[0])
    signalvar_static = emmap1.signalvar_lst[0]
    # smoothening signal
    fobj.write("\nStatic map signal \n")
    signalvar_static = get_extended_signal(
        res_arr, mapaverage.set_array(signalvar_static), hffsc_lst[0], fobj=fobj
    )
    signalvar_lst.append(mapaverage.set_array(signalvar_static))
    i = -1
    for fo, t, rotmat in zip(fo_lst[1:], t_lst, r_lst):
        i = i + 1
        st, _, _, _ = fcodes_fast.get_st(nx, ny, nz, t)
        frt_full = utils.get_FRS(cell, rotmat, fo * st)[:, :, :, 0]
        data2write = np.real(np.fft.ifftshift(np.fft.ifftn(np.fft.ifftshift(frt_full))))
        data2write = shift(data2write, np.subtract(com_lst[0], emmap1.box_centr))
        core.iotools.write_mrc(
            data2write, "fitted_map_" + str(i + 1) + "_.mrc", cell, origin
        )
        # estimating covaraince between current map vs. static map
        _, f1f2_fsc = fcodes_fast.calc_covar_and_fsc_betwn_anytwomaps(
            fo_lst[0], frt_full, bin_idx, nbin, 0, nx, ny, nz
        )
        # find the relative B-factor between static and frt_full
        resol = mapaverage.get_resol(f1f2_fsc, res_arr)
        frt_hf1 = utils.get_FRS(cell, rotmat, emmap1.fhf_lst[2 * i + 2] * st)[
            :, :, :, 0
        ]
        frt_hf2 = utils.get_FRS(cell, rotmat, emmap1.fhf_lst[2 * i + 3] * st)[
            :, :, :, 0
        ]
        s_grid = (emmap1.s_grid ** 2) / 4.0
        scale1, bfac_relative1 = mapaverage.estimate_scale_and_relativebfac(
            [frt_lst[0], frt_hf1], cell, resol
        )
        scale2, bfac_relative2 = mapaverage.estimate_scale_and_relativebfac(
            [frt_lst[0], frt_hf2], cell, resol
        )
        frt_hf1_scaled = scale1 * frt_hf1 * np.exp(bfac_relative1 * s_grid)
        frt_hf2_scaled = scale2 * frt_hf2 * np.exp(bfac_relative2 * s_grid)
        # re-estimate sigvar and totvar
        (
            bin_fsc,
            _,
            signalvar,
            totalvar,
            frt_full_scaled,
            _,
        ) = core.fsc.halfmaps_fsc_variance(
            frt_hf1_scaled, frt_hf2_scaled, bin_idx, nbin
        )
        maps2average_lst.append(frt_full_scaled)
        # smoothen signal
        bin_fsc_full = 2 * bin_fsc / (1 + bin_fsc)
        fobj.write("\nMoving map signal \n")
        signalvar = get_extended_signal(
            res_arr, mapaverage.set_array(signalvar), bin_fsc_full, fobj=fobj
        )
        signalvar_lst.append(mapaverage.set_array(signalvar))

        f1f2_covar, f1f2_fsc = fcodes_fast.calc_covar_and_fsc_betwn_anytwomaps(
            fo_lst[0], frt_full_scaled, bin_idx, nbin, 0, nx, ny, nz
        )
        # smoothen signal
        fobj.write("\nCovariance between static and moving maps \n")
        f1f2_covar = get_extended_signal(
            res_arr, mapaverage.set_array(f1f2_covar), f1f2_fsc, fobj=fobj
        )
        # mask f1f2_covar so that anything beyond zero get zeros
        cov_lst.append(mapaverage.set_array(f1f2_covar))
        # signalvar_lst.append(mapaverage.set_array(signalvar))
        totalvar_lst.append(mapaverage.set_array(totalvar))

    var_cov_lst = signalvar_lst + cov_lst
    core.plotter.plot_nlines_log(res_arr, var_cov_lst)
    # map averaging with un-normalised maps
    averagemaps = utils.calc_averagemaps_simple(
        maps2average_lst,
        cell,
        nbin,
        signalvar_lst,
        totalvar_lst,
        cov_lst,
        bin_idx,
        emmap1.s_grid,
        res_arr,
        Bf_arr,
        emmap1.com_lst,
        emmap1.box_centr,
        origin,
    )

    utils.output_maps(
        averagemaps,
        emmap1.com_lst,
        t_lst,
        r_lst,
        emmap1.box_centr,
        cell,
        origin,
        Bf_arr,
    )


def output_rotated_maps_using_fo_allmaps(emmap1, r_lst, t_lst, Bf_arr=None, fobj=None):
    from scipy.ndimage.interpolation import shift

    if Bf_arr is None:
        Bf_arr = [0.0]
    com_lst = emmap1.com_lst
    fo_lst = emmap1.fo_lst
    nmaps = len(fo_lst)
    signalvar_lst = []
    totalvar_lst = []
    totalvar_lst.append(emmap1.totalvar_lst[0])
    hffsc_lst = []
    hffsc_lst.append(mapaverage.set_array(emmap1.hffsc_lst[0]))
    bin_idx = emmap1.bin_idx
    nbin = emmap1.nbin
    res_arr = emmap1.res_arr
    cell = emmap1.map_unit_cell
    origin = emmap1.map_origin
    nx, ny, nz = fo_lst[0].shape
    frt_lst = []
    frt_lst.append(fo_lst[0])
    ert_lst = []
    ert_lst.append(emmap1.eo_lst[0])
    cov_lst = []
    var_covar_mat = np.zeros((nmaps, nmaps, nbin), dtype="float")
    # static map
    data2write = np.real(np.fft.ifftshift(np.fft.ifftn(np.fft.ifftshift(fo_lst[0]))))
    data2write = shift(data2write, np.subtract(com_lst[0], emmap1.box_centr))
    core.iotools.write_mrc(data2write, "static_map.mrc", cell, origin)
    del data2write
    fobj.write("\n *** Map averaging using un-normalized structure factors *** \n")
    maps2average_lst = []
    maps2average_lst.append(fo_lst[0])
    signalvar_static = emmap1.signalvar_lst[0]
    signalvar_lst.append(mapaverage.set_array(signalvar_static))
    var_covar_mat[0, 0, :] = mapaverage.set_array(signalvar_static)
    s_grid = (emmap1.s_grid ** 2) / 4.0
    i = 0
    for fo, t, rotmat in zip(fo_lst[1:], t_lst, r_lst):
        i = i + 1
        st, _, _, _ = fcodes_fast.get_st(nx, ny, nz, t)
        frt_full = utils.get_FRS(cell, rotmat, fo * st)[:, :, :, 0]
        data2write = np.real(np.fft.ifftshift(np.fft.ifftn(np.fft.ifftshift(frt_full))))
        data2write = shift(data2write, np.subtract(com_lst[0], emmap1.box_centr))
        core.iotools.write_mrc(
            data2write, "fitted_map_" + str(i) + "_.mrc", cell, origin
        )
        # estimating covaraince between current map vs. static map
        _, f1f2_fsc = fcodes_fast.calc_covar_and_fsc_betwn_anytwomaps(
            fo_lst[0], frt_full, bin_idx, nbin, 0, nx, ny, nz
        )
        # get resolution from fsc
        resol = mapaverage.get_resol(f1f2_fsc, res_arr)
        # apply transformation on half maps of moving map
        frt_hf1 = utils.get_FRS(cell, rotmat, emmap1.fhf_lst[2 * i] * st)[:, :, :, 0]
        frt_hf2 = utils.get_FRS(cell, rotmat, emmap1.fhf_lst[2 * i + 1] * st)[
            :, :, :, 0
        ]
        frt_lst.append((frt_hf1 + frt_hf2) / 2.0)
        # find the relative B-factor between static and frt_full
        scale1, bfac_relative1 = mapaverage.estimate_scale_and_relativebfac(
            [frt_lst[0], frt_hf1], cell, resol
        )
        scale2, bfac_relative2 = mapaverage.estimate_scale_and_relativebfac(
            [frt_lst[0], frt_hf2], cell, resol
        )
        frt_hf1_scaled = scale1 * frt_hf1 * np.exp(bfac_relative1 * s_grid)
        frt_hf2_scaled = scale2 * frt_hf2 * np.exp(bfac_relative2 * s_grid)
        # re-estimate sigvar and totvar
        (
            _,
            _,
            signalvar,
            totalvar,
            frt_full_scaled,
            ert_full_scaled,
        ) = core.fsc.halfmaps_fsc_variance(
            frt_hf1_scaled, frt_hf2_scaled, bin_idx, nbin
        )
        maps2average_lst.append(frt_full_scaled)
        ert_lst.append(ert_full_scaled)
        signalvar_lst.append(mapaverage.set_array(signalvar))
        var_covar_mat[i, i, :] = mapaverage.set_array(signalvar)

        f1f2_covar, f1f2_fsc = fcodes_fast.calc_covar_and_fsc_betwn_anytwomaps(
            fo_lst[0], frt_full_scaled, bin_idx, nbin, 0, nx, ny, nz
        )
        # mask f1f2_covar so that anything beyond zero get zeros
        cov_lst.append(mapaverage.set_array(f1f2_covar))
        var_covar_mat[0, i, :] = mapaverage.set_array(f1f2_covar)
        var_covar_mat[i, 0, :] = mapaverage.set_array(f1f2_covar)
        # signalvar_lst.append(mapaverage.set_array(signalvar))
        totalvar_lst.append(mapaverage.set_array(totalvar))

    # calculating covariances between fullmaps
    fobj.write("\nCovariance between moving maps after fitting! \n")
    t_init = [0.0, 0.0, 0.0]
    theta_init = [(1, 0, 0), 0.0]
    ncycles = 10
    smax = 6
    for i in range(len(frt_lst)):
        for j in range(len(frt_lst)):
            if i == j:
                continue
            if i == 0 or j == 0:
                continue
            if j > i:
                frt_i = frt_lst[i]
                frt_j = frt_lst[j]
                # Can you do fitting between i and j here before calc. covariance
                emmap1.eo_lst = [ert_lst[i], ert_lst[j]]
                q = core.quaternions.get_quaternion(theta_init)
                q = q / np.sqrt(np.dot(q, q))
                rotmat = core.quaternions.get_RM(q)
                if len(emmap1.res_arr) < 50:
                    slf = len(emmap1.res_arr)
                else:
                    slf = 50
                rotmat, trans = run_fit.run_fit(
                    emmap1=emmap1,
                    smax=smax,
                    rotmat=rotmat,
                    t=t_init,
                    slf=slf,
                    ncycles=ncycles,
                    ifit=0,
                    fobj=fobj,
                    interp="linear",
                )
                st, _, _, _ = fcodes_fast.get_st(nx, ny, nz, trans)
                frt_j_new = utils.get_FRS(cell, rotmat, frt_j * st)[:, :, :, 0]
                # scale 1st map to 2nd map
                _, f1f2_fsc = fcodes_fast.calc_covar_and_fsc_betwn_anytwomaps(
                    frt_i, frt_j_new, bin_idx, nbin, 0, nx, ny, nz
                )
                # get resolution from fsc
                resol = mapaverage.get_resol(f1f2_fsc, res_arr)
                scale1, bfac_relative1 = mapaverage.estimate_scale_and_relativebfac(
                    [frt_i, frt_j_new], cell, resol
                )
                frt_j_scaled = scale1 * frt_j_new * np.exp(bfac_relative1 * s_grid)
                # calc covariance
                f1f2_covar, f1f2_fsc = fcodes_fast.calc_covar_and_fsc_betwn_anytwomaps(
                    frt_i, frt_j_scaled, bin_idx, nbin, 0, nx, ny, nz
                )
                # smoothen signal
                fobj.write("\nCovariance between map_i and map_j maps \n")
                f1f2_covar = get_extended_signal(
                    res_arr, mapaverage.set_array(f1f2_covar), f1f2_fsc, fobj=fobj
                )
                # mask f1f2_covar so that anything beyond zero get zeros
                cov_lst.append(mapaverage.set_array(f1f2_covar))
                var_covar_mat[i, j, :] = mapaverage.set_array(f1f2_covar)

    # map averaging with un-normalised maps
    averagemaps = utils.calc_averagemaps_simple_allmaps(
        maps2average_lst,
        cell,
        nbin,
        var_covar_mat,
        totalvar_lst,
        bin_idx,
        emmap1.s_grid,
        res_arr,
        Bf_arr,
    )

    # calc. fsc between static map and averaged maps
    utils.output_maps(
        averagemaps,
        emmap1.com_lst,
        t_lst,
        r_lst,
        emmap1.box_centr,
        cell,
        origin,
        Bf_arr,
    )


def avgmaps_using_fo_nofit(emmap1, Bf_arr=None, fobj=None):
    if Bf_arr is None:
        Bf_arr = [0.0]
    fo_lst = emmap1.fo_lst
    nmaps = len(fo_lst)
    signalvar_lst = []
    totalvar_lst = []
    totalvar_lst.append(mapaverage.set_array(emmap1.totalvar_lst[0]))
    hffsc_lst = []
    hffsc_lst.append(mapaverage.set_array(emmap1.hffsc_lst[0]))
    bin_idx = emmap1.bin_idx
    nbin = emmap1.nbin
    res_arr = emmap1.res_arr
    cell = emmap1.map_unit_cell
    origin = emmap1.map_origin
    nx, ny, nz = fo_lst[0].shape
    frt_lst = []
    frt_lst.append(fo_lst[0])
    ert_lst = []
    ert_lst.append(emmap1.eo_lst[0])
    var_covar_mat = np.zeros((nmaps, nmaps, nbin), dtype="float")
    # static map
    data2write = np.real(np.fft.ifftshift(np.fft.ifftn(np.fft.ifftshift(fo_lst[0]))))
    core.iotools.write_mrc(data2write, "static_map.mrc", cell, origin)
    del data2write
    fobj.write("\n *** Map averaging using un-normalized structure factors *** \n")
    maps2average_lst = []
    maps2average_lst.append(emmap1.fo_lst[0])
    signalvar_static = emmap1.signalvar_lst[0]
    signalvar_lst.append(mapaverage.set_array(signalvar_static))
    var_covar_mat[0, 0, :] = mapaverage.set_array(signalvar_static)
    i = 0
    for i in range(fo_lst[1:]):
        i = i + 1
        # scaling using power spectrum
        from emda.ext.scale_maps import scale_twomaps_by_power, transfer_power

        scale_grid = transfer_power(
            bin_idx=bin_idx,
            res_arr=res_arr,
            scale=scale_twomaps_by_power(
                f1=emmap1.fo_lst[0],
                f2=emmap1.fo_lst[i],
                bin_idx=bin_idx,
                uc=cell,
                res_arr=res_arr,
            ),
        )
        fhf1 = emmap1.fhf_lst[2 * i] * scale_grid
        fhf2 = emmap1.fhf_lst[2 * i + 1] * scale_grid
        (_, _, signalvar, totalvar, fo_scaled, _,) = core.fsc.halfmaps_fsc_variance(
            fhf1, fhf2, bin_idx, nbin
        )
        var_covar_mat[i, i, :] = mapaverage.set_array(signalvar)
        totalvar_lst.append(mapaverage.set_array(totalvar))
        maps2average_lst.append(fo_scaled)
        f1f2_covar, _, _ = fcodes_fast.calc_covar_and_fsc_betwn_anytwomaps(
            fo_lst[0], fo_scaled, bin_idx, nbin, 0, nx, ny, nz
        )
        var_covar_mat[0, i, :] = mapaverage.set_array(f1f2_covar)
        var_covar_mat[i, 0, :] = mapaverage.set_array(f1f2_covar)

    # map averaging with un-normalised maps
    averagemaps = utils.calc_averagemaps_simple_allmaps(
        maps2average_lst,
        cell,
        nbin,
        var_covar_mat,
        totalvar_lst,
        bin_idx,
        emmap1.s_grid,
        res_arr,
        Bf_arr,
    )
    return averagemaps


def get_extended_signal(res_arr, signal, bin_fsc, fobj=None):
    from scipy import stats
    from emda.ext.mapfit import rdp_algo

    fobj.write("\nSignal extended using B-factor \n")
    bin_fsc_masked = bin_fsc[bin_fsc > 0.0]
    res_arr_masked = res_arr[bin_fsc > 0.0]
    dist2 = np.sqrt((bin_fsc_masked - 0.3) ** 2)
    uplim = np.argmin(dist2)
    # Need to determine the low and high limits
    x = 1.0 / res_arr_masked[:uplim]
    y = signal[:uplim]
    xc, _ = rdp_algo.run_rdp(x, np.log(y), epsilon=0.01)
    upl = np.argmin(np.sqrt((x - xc[-1]) ** 2))
    lwl = np.argmin(np.sqrt((x - xc[-2]) ** 2))
    s1 = x[lwl:upl]
    sig1 = y[lwl:upl]
    fobj.write("low resol. limit: " + str(res_arr[signal > 0.0][lwl]) + " (A) \n")
    fobj.write("high resol. limit: " + str(res_arr[signal > 0.0][upl]) + " (A) \n")
    slope, intercept, _, _, _ = stats.linregress((s1 * s1) / 4.0, np.log(sig1))
    bfac = slope
    print("B factor = ", bfac, "Intercept = ", intercept)
    fobj.write("B-factor: " + str(bfac) + " \n")
    fobj.write("Intercept: " + str(intercept) + " \n")
    s = 1.0 / res_arr
    signal_pred = np.exp((bfac / 4) * s ** 2 + intercept)
    new_signal = signal
    new_signal[upl:] = signal_pred[upl:]
    return new_signal
