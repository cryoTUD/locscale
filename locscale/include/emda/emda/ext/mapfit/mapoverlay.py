"""
Author: "Rangana Warshamanage, Garib N. Murshudov"
MRC Laboratory of Molecular Biology
    
This software is released under the
Mozilla Public License, version 2.0; see LICENSE.
"""

from __future__ import absolute_import, division, print_function, unicode_literals
from timeit import default_timer as timer
import numpy as np
import fcodes_fast
from emda import core
from emda.ext.mapfit import map_class
from numpy.fft import fftn, ifftn, fftshift, ifftshift
from emda.ext.overlay import output_rotated_maps as outputmaps
from emda.ext.overlay import output_rotated_models


def output_rotated_maps(emmap1, r_lst, t_lst, Bf_arr=None):
    from emda.ext.mapfit import utils
    from scipy.ndimage.interpolation import shift

    if Bf_arr is None:
        Bf_arr = [0.0]
    #emmap1.com = False
    if emmap1.com:
        com = emmap1.com1
    fo_lst = emmap1.fo_lst
    bin_idx = emmap1.bin_idx
    nbin = emmap1.nbin
    res_arr = emmap1.res_arr
    cell = emmap1.map_unit_cell
    origin = emmap1.map_origin
    nx, ny, nz = fo_lst[0].shape
    frt_lst = []
    frt_lst.append(fo_lst[0])
    cov_lst = []
    fsc12_lst = []
    fsc12_lst_unaligned = []
    imap_f = 0
    # static map
    data2write = np.real(ifftshift(ifftn(ifftshift(fo_lst[0]))))
    if emmap1.com:
        data2write = shift(data2write, np.subtract(com, emmap1.box_centr))
    core.iotools.write_mrc(data2write, "static_map.mrc", cell, origin)
    del data2write
    for fo, t, rotmat in zip(fo_lst[1:], t_lst, r_lst):
        f1f2_fsc_unaligned = core.fsc.anytwomaps_fsc_covariance(
            fo_lst[0], fo, bin_idx, nbin
        )[0]
        fsc12_lst_unaligned.append(f1f2_fsc_unaligned)
        imap_f = imap_f + 1
        st, _, _, _ = fcodes_fast.get_st(nx, ny, nz, t)
        frt_full = utils.get_FRS(rotmat, fo * st, interp="cubic")[:, :, :, 0]
        frt_lst.append(frt_full)
        data2write = np.real(ifftshift(ifftn(ifftshift(frt_full))))
        if emmap1.com:
            data2write = shift(data2write, np.subtract(com, emmap1.box_centr))
        core.iotools.write_mrc(
            data2write,
            "{0}_{1}.{2}".format("fitted_map", str(imap_f), "mrc"),
            cell,
            origin,
        )
        # estimating covaraince between current map vs. static map
        bin_stats = core.fsc.anytwomaps_fsc_covariance(
            fo_lst[0], frt_full, bin_idx, nbin
        )
        f1f2_fsc, f1f2_covar = bin_stats[0], bin_stats[1]
        cov_lst.append(f1f2_covar)
        fsc12_lst.append(f1f2_fsc)
        core.plotter.plot_nlines(
            res_arr,
            [f1f2_fsc_unaligned, f1f2_fsc],
            "{0}_{1}.{2}".format("fsc", str(imap_f), "eps"),
            ["FSC before", "FSC after"],
        )


def main(
    maplist,
    ncycles,
    t_init,
    theta_init,
    smax,
    masklist,
    fobj,
    interp,
    modelres=5.0,
    halfmaps=False,
    dfs_interp=False,
    usemodel=False,
    usecom=False,
    fitres=None,
):
    from emda.ext.mapfit import utils, run_fit, interp_derivatives

    if len(maplist) < 2:
        print(" At least 2 maps required!")
        exit()
    try:
        if halfmaps:
            print("Map overlay using halfmaps")
            emmap1 = map_class.Overlay(maplist, masklist)
        else:
            print("Map overlay not using halfmaps")
            emmap1 = map_class.EmmapOverlay(map_list=maplist, mask_list=masklist, modelres=modelres, com=usecom)
        fobj.write("Map overlay\n")
    except NameError:
        if halfmaps:
            print("Map overlay using halfmaps")
            emmap1 = map_class.Overlay(map_list=maplist)
        else:
            print("Map overlay not using halfmaps")
            emmap1 = map_class.EmmapOverlay(map_list=maplist, modelres=modelres, com=usecom)
        fobj.write("Map overlay\n")
    if usemodel:
        emmap1.load_models()
    else:
        emmap1.load_maps(fobj)
    emmap1.calc_fsc_from_maps(fobj)
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
    t = [itm / emmap1.pixsize[i] for i, itm in enumerate(t_init)]
    rotmat_lst = []
    transl_lst = []
    # resolution estimate for line-fit
    dist = np.sqrt((emmap1.res_arr - smax) ** 2)
    slf = np.argmin(dist) + 1
    if slf % 2 != 0:
        slf = slf - 1
    slf = min([len(dist), slf])

    for ifit in range(1, len(emmap1.eo_lst)):
        fobj.write("Fitting set#: " + str(ifit) + " \n")
        start_fit = timer()
        if dfs_interp:
            cell = emmap1.map_unit_cell
            xyz = utils.create_xyz_grid(cell, emmap1.map_dim)
            vol = cell[0] * cell[1] * cell[2]
            dfs_full = interp_derivatives.dfs_fullmap(
                ifftshift(np.real(ifftn(ifftshift(emmap1.eo_lst[ifit])))),
                xyz,
                vol,
            )
        else:
            dfs_full = None
        rotmat, translation = run_fit.run_fit(
            emmap1=emmap1,
            smax=smax,
            rotmat=rotmat,
            t=t,
            slf=slf,
            ncycles=ncycles,
            ifit=ifit,
            fobj=fobj,
            interp=interp,
            dfs_full=dfs_full,
            fitres=fitres,
        )
        rotmat_lst.append(rotmat)
        transl_lst.append(translation)
        end_fit = timer()
        fobj.write("Final Translation: " + str(translation) + " \n")
        fobj.write("Final Rotation matrix: " + str(rotmat) + " \n")
        print("time for fitting: ", end_fit - start_fit)

    #output_rotated_maps(emmap1, rotmat_lst, transl_lst)
    outputmaps(emmap1, r_lst=rotmat_lst, t_lst=transl_lst)
    output_rotated_models(emmap1, maplist=maplist, r_lst=rotmat_lst, t_lst=transl_lst)
