# Difference map using half data
# 1. read half maps for each
# 2. calculate bestmap for each
# 3. compute difference map from best maps
import numpy as np
import fcodes_fast
from emda.ext import bestmap
from emda.core import iotools, restools, quaternions, fsc, plotter
from numpy.fft import fftn, ifftn, fftshift, ifftshift
from emda.ext.overlay import EmmapOverlay, run_fit
from emda.ext.mapfit import utils
import emda.emda_methods as em


mode = 1 # bestmap in resolution bins
#mode = 2 # bestmap in 3d
maplist = [
        "/Users/ranganaw/MRC/REFMAC/Vinoth/nat_half1_class001_unfil.map",
        "/Users/ranganaw/MRC/REFMAC/Vinoth/nat_half2_class001_unfil.map",
        "/Users/ranganaw/MRC/REFMAC/Vinoth/lig_half1_class001_unfil.map",
        "/Users/ranganaw/MRC/REFMAC/Vinoth/lig_half2_class001_unfil.map",
]

def get_bestmap(arr1, arr2, bin_idx, nbin, mode):
    f1 = fftshift(fftn(arr1))
    f2 = fftshift(fftn(arr2))
    f_map = bestmap.bestmap(f1=f1, f2=f2, bin_idx=bin_idx, nbin=nbin, mode=mode)
    return f_map

def fitmaps(emmap1, smax=6, ncy=5, fitres=None):
    rotmat_init = np.identity(3)
    t_init=[0.0, 0.0, 0.0]
    t = [itm / emmap1.pixsize for itm in t_init]
    # resolution estimate for line-fit
    dist = np.sqrt((emmap1.res_arr - smax) ** 2)
    slf = np.argmin(dist) + 1
    if slf % 2 != 0:
        slf = slf - 1
    slf = min([len(dist), slf])
    rotmat_list = []
    trans_list = []
    for ifit in range(1, len(emmap1.eo_lst)):
        t, q_final = run_fit(
            emmap1=emmap1,
            smax=smax,
            rotmat=rotmat_init,
            t=t,
            slf=slf,
            ncycles=ncy,
            ifit=ifit,
            interp="linear",
            fitres=fitres,
        )
        rotmat = quaternions.get_RM(q_final)
        rotmat_list.append(rotmat)
        trans_list.append(t)
    # output maps
    return emmap1, rotmat_list, trans_list  

def output_rotated_maps(emmap1, r_lst, t_lst, bin_diffmap):
    fo_lst = emmap1.fo_lst
    bin_idx = emmap1.bin_idx
    res_arr = emmap1.res_arr
    nbin = emmap1.nbin
    f_static = fo_lst[0]
    nx, ny, nz = f_static.shape
    fout1 = restools.cut_resolution(f_static, bin_idx, res_arr, bin_diffmap)
    i = 0
    for fo, t, rotmat in zip(fo_lst[1:], t_lst, r_lst):
        i += 1
        f1f2_fsc_unaligned = fsc.anytwomaps_fsc_covariance(
            f_static, fo, bin_idx, nbin
        )[0]
        f1f2_fsc_unaligned = np.nan_to_num(f1f2_fsc_unaligned, copy=False, nan=0.0) # for aesthetics
        st, _, _, _ = fcodes_fast.get_st(nx, ny, nz, t)
        frt = utils.get_FRS(rotmat, fo * st, interp="cubic")[:, :, :, 0]
        # estimating covaraince between current map vs. static map
        f1f2_fsc = fsc.anytwomaps_fsc_covariance(f_static, frt, bin_idx, nbin)[0]
        f1f2_fsc = np.nan_to_num(f1f2_fsc, copy=False, nan=0.0) # for aesthetics
        plotter.plot_nlines(
            emmap1.res_arr,
            [f1f2_fsc_unaligned[: emmap1.nbin], f1f2_fsc[: emmap1.nbin]],
            "{0}_{1}.{2}".format("fsc", str(i), "eps"),
            ["FSC before", "FSC after"],
        )
        fout2 = restools.cut_resolution(frt, bin_idx, res_arr, bin_diffmap)

    diffm1m2 = fout1 - fout2
    diffm2m1 = fout2 - fout1
    diffmap = np.stack((diffm1m2, diffm2m1, fout1, fout2), axis = -1)
    return diffmap


def check_sampling(arr, tpix, cpix, target_dim):
    resampled_arr = iotools.resample2staticmap(
        curnt_pix=cpix,
        targt_pix=tpix,
        targt_dim=target_dim,
        arr=arr,
    )
    return resampled_arr


def main(maplist, results, fit=True, resol=3, ncy=5, usecom=False, fitres=None, masklist=None):
    bestmap_list = []
    msk_list = []
    nmaps = len(maplist)
    if nmaps == 3:
        print("TO DO.")
        """ uc, arr1, origin = iotools.read_map(maplist[0])
        uc, arr2, origin = iotools.read_map(maplist[1])
        target_uc = uc
        target_dim = arr1.shape
        tpix = target_uc[0] / target_dim[0]
        nbin, res_arr, bin_idx = restools.get_resolution_array(uc, arr1)
        bestmap_list.append(get_bestmap(arr1, arr2, bin_idx, nbin, mode))
        if maplist[2].endswith(((".pdb", ".cif", ".ent"))):
            # calculate map from model
            modelmap = em.model2map(
                modelxyz=maplist[2],
                dim=arr1.shape,
                resol=resol,
                cell=uc,
                maporigin=origin,
            )
            bestmap_list.append(fftshift(fftn(modelmap))) """
    elif nmaps == 4:
        if masklist is not None and len(maplist) // 2 != len(masklist):
            print("Number of masks should be equal \
                to half of the number of maps")
            print("Masks will not be used in the calculation")
            masklist = None
        for i in range(0, nmaps, 2):
            uc, arr1, origin = iotools.read_map(maplist[i])
            uc, arr2, origin = iotools.read_map(maplist[i+1])
            if masklist is not None:
                _, msk, _ = iotools.read_map(masklist[i//2])
                try:
                    assert arr1.shape == msk.shape
                    msk_list.append(msk)
                except AssertionError as error:
                    print(error)
            if i == 0:
                target_uc = uc
                map_origin = origin
                target_dim = arr1.shape
                tpix = target_uc[0] / target_dim[0]
                nbin, res_arr, bin_idx = restools.get_resolution_array(uc, arr1)
            else:
                cpix = uc[0] / arr1.shape[0]
            bestmap_list.append(get_bestmap(arr1, arr2, bin_idx, nbin, mode))

    # check for sampling in best maps
    arr_to_sample = np.real(ifftn(ifftshift(bestmap_list[1])))
    resampled_arr = check_sampling(arr_to_sample, tpix, cpix, target_dim)
    if len(msk_list) == 2:
        resampled_msk = check_sampling(msk_list[1], tpix, cpix, target_dim)
        msk_list[1] = resampled_msk
        results.masklist = msk_list
    uc = target_uc
    results.cell = uc
    results.origin = map_origin
    # maps for fitting
    fmaps_to_fit = bestmap_list
    if len(msk_list) == 2:
        static_bestmap = np.real(ifftn(ifftshift(bestmap_list[0])))
        fmaps_to_fit[0] = fftshift(fftn(static_bestmap * msk_list[0]))
        fmaps_to_fit[1] = fftshift(fftn(resampled_arr * msk_list[1]))
    else:
        fmaps_to_fit[1] = fftshift(fftn(resampled_arr))

    bin_diffmap = np.argmin(np.sqrt((res_arr - resol) ** 2)) + 1

    if fit:
        # prepare bestmaps for fitting
        print("*** Difference map with fit optimization ***")
        maps = []
        for j, imap in enumerate(fmaps_to_fit):
            if j == 0:
                nx, ny, nz = imap.shape
                t = np.array([0.5, 0.5, 0.5], dtype='float')
                st, _, _, _ = fcodes_fast.get_st(nx, ny, nz, t)
            maps.append(imap * st)
        # TO DO: INCLUDE USECOM FOR FIT
        emmap1 = EmmapOverlay(map_list=maplist)
        emmap1.map_dim = maps[0].shape
        emmap1.map_unit_cell = uc
        emmap1.pixsize = uc[0] / emmap1.map_dim[0]
        emmap1.map_origin = map_origin
        emmap1.nbin = nbin
        emmap1.bin_idx = bin_idx
        emmap1.res_arr = res_arr
        emmap1.fo_lst = maps
        emmap1.eo_lst = maps
        emmap1, rotmat_list, trans_list = fitmaps(emmap1, ncy=ncy, fitres=fitres)
        diffmap = output_rotated_maps(emmap1, rotmat_list, trans_list, bin_diffmap)
    else:
        # difference map without fit
        fout1 = restools.cut_resolution(fmaps_to_fit[0], bin_idx, res_arr, bin_diffmap)
        fout2 = restools.cut_resolution(fmaps_to_fit[1], bin_idx, res_arr, bin_diffmap)
        diffm1m2 = fout1 - fout2
        diffm2m1 = fout2 - fout1
        diffmap = np.stack((diffm1m2, diffm2m1, fout1, fout2), axis = -1)
    results.diffmap = diffmap
    return results


if __name__ == "__main__":
    results = main(maplist, results)