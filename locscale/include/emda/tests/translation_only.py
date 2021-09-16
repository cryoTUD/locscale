"""
Author: "Rangana Warshamanage, Garib N. Murshudov"
MRC Laboratory of Molecular Biology
    
This software is released under the
Mozilla Public License, version 2.0; see LICENSE.
"""

from __future__ import absolute_import, division, print_function, unicode_literals

from timeit import default_timer as timer
import math
import numpy as np
import fcodes_fast
import emda.emda_methods as em
from numpy.fft import fftn, ifftn, fftshift, ifftshift
from emda import core
from emda.ext.mapfit import utils
from emda.ext.overlay import cut_resolution_for_linefit



np.set_printoptions(suppress=True)  # Suppress insignificant values for clarity

timeit = False


class EmmapOverlay:
    def __init__(self, map_list, mask_list=None):
        self.map_list = map_list
        self.mask_list = mask_list
        self.map_unit_cell = None
        self.map_origin = None
        self.map_dim = None
        self.pixsize = None
        self.arr_lst = []
        self.ceo_lst = None
        self.cfo_lst = None
        self.cbin_idx = None
        self.cdim = None
        self.cbin = None
        self.box_centr = None
        self.fhf_lst = None
        self.nbin = None
        self.res_arr = None
        self.bin_idx = None
        self.fo_lst = None
        self.eo_lst = None
        self.totalvar_lst = None

    def load_maps(self):
        from scipy import ndimage
        from scipy.ndimage.interpolation import shift

        cmask = False
        fhf_lst = []
        full_fhf_lst = []
        if self.mask_list is not None:
            if len(self.map_list) != len(self.mask_list):
                raise SystemExit("map_list and mask_list must have the same size!")
            for i in range(len(self.mask_list)):
                if i == 0:
                    _, mask, _ = em.get_data(self.mask_list[i])
                    mask = utils.set_dim_even(mask)
                    uc, arr, origin = em.get_data(self.map_list[i])
                    arr = utils.set_dim_even(arr)
                    try:
                        assert arr.shape == mask.shape
                    except AssertionError:
                        raise SystemExit("Map and Mask Dimension mismatched!")
                    arr = arr * mask
                    nx, ny, nz = arr.shape
                    map_origin = origin
                    uc_target = uc
                    target_dim = arr.shape
                    target_pix_size = uc_target[0] / target_dim[0]
                    if cmask:
                        corner_mask = utils.remove_unwanted_corners(uc, target_dim)
                    else:
                        corner_mask = 1.0
                    fhf_lst.append(fftshift(fftn(fftshift(arr * corner_mask))))
                else:
                    uc, arr, origin = em.get_data(
                        self.map_list[i],
                        dim=target_dim,
                        uc=uc_target,
                        maporigin=map_origin,
                    )
                    arr = utils.set_dim_even(arr)
                    print("origin: ", origin)
                    _, mask, _ = em.get_data(self.mask_list[i])
                    mask = utils.set_dim_even(mask)
                    try:
                        assert arr.shape == mask.shape
                    except AssertionError:
                        raise SystemExit("Map and Mask Dimension mismatched!")
                    arr = arr * mask
                    curnt_pix_size = uc[0] / arr.shape[0]
                    arr = core.iotools.resample2staticmap(
                        curnt_pix=curnt_pix_size,
                        targt_pix=target_pix_size,
                        targt_dim=target_dim,
                        arr=arr,
                    )
                    self.arr_lst.append(arr * corner_mask)
                    fhf_lst.append(fftshift(fftn(fftshift(arr * corner_mask))))
            self.pixsize = target_pix_size
            self.map_origin = map_origin
            self.map_unit_cell = uc_target
            self.map_dim = target_dim
            self.fhf_lst = fhf_lst
        if self.mask_list is None:
            for i in range(len(self.map_list)):
                if i == 0:
                    uc, arr, origin = em.get_data(self.map_list[i])
                    arr = utils.set_dim_even(arr)
                    print("origin: ", origin)
                    nx, ny, nz = arr.shape
                    map_origin = origin
                    uc_target = uc
                    target_dim = arr.shape
                    target_pix_size = []
                    for i in range(3):
                        target_pix_size.append(uc_target[i] / target_dim[i])
                    self.arr_lst.append(arr)
                    fhf_lst.append(fftshift(fftn(fftshift(arr))))
                else:
                    uc, arr, origin = em.get_data(
                        self.map_list[i],
                        dim=target_dim,
                        uc=uc_target,
                        maporigin=map_origin,
                    )
                    arr = utils.set_dim_even(arr)
                    em.write_mrc(arr, 'modelmap'+str(i)+'.mrc', uc, origin)
                    print("origin: ", origin)
                    curnt_pix_size = []
                    for i in range(3):
                        curnt_pix_size.append(uc[i] / arr.shape[i])
                    arr = core.iotools.resample2staticmap(
                        curnt_pix=curnt_pix_size,
                        targt_pix=target_pix_size,
                        targt_dim=target_dim,
                        arr=arr,
                    )
                    self.arr_lst.append(arr)
                    fhf_lst.append(fftshift(fftn(fftshift(arr))))
            self.pixsize = target_pix_size
            self.map_origin = map_origin
            self.map_unit_cell = uc_target
            self.map_dim = target_dim
            self.fhf_lst = fhf_lst

    def calc_fsc_from_maps(self):
        # function for only two maps fitting
        nmaps = len(self.fhf_lst)
        fFo_lst = []
        fEo_lst = []
        fBTV_lst = []

        self.nbin, self.res_arr, self.bin_idx = core.restools.get_resolution_array(
            self.map_unit_cell, self.fhf_lst[0]
        )
        #
        for i in range(nmaps):
            _, _, _, totalvar, fo, eo = core.fsc.halfmaps_fsc_variance(
                self.fhf_lst[i], self.fhf_lst[i], self.bin_idx, self.nbin
            )
            fFo_lst.append(fo)
            fEo_lst.append(eo)
            fBTV_lst.append(totalvar)
        #
        self.fo_lst = fFo_lst
        self.eo_lst = fEo_lst
        self.totalvar_lst = fBTV_lst

def derivatives_translation(e0, e1, wgrid, w2grid, sv):
    PI = np.pi
    tp2 = (2.0 * PI)**2
    tpi = (2.0 * PI * 1j)
    # translation derivatives
    df = np.zeros(3, dtype='float')
    ddf = np.zeros((3,3), dtype='float')
    for i in range(3):
        df[i] = np.real(np.sum(wgrid * e0 * np.conjugate(e1 * tpi * sv[i,:,:,:])))
        for j in range(3):
            if(i==0 or (i>0 and j>=i)):
                ddf[i,j] = -tp2 * np.sum(w2grid * sv[i,:,:,:] * sv[j,:,:,:])
            else:
                ddf[i,j] = ddf[j,i]
    ddf_inv = np.linalg.pinv(ddf)
    step = ddf_inv.dot(-df)
    return step

def trans2angstroms(t, pixsize, dim):
    t_angstrom = np.asarray(t, 'float') * np.asarray(pixsize, 'float') * np.asarray(dim, 'int')
    return np.sqrt(np.dot(t_angstrom, t_angstrom))

class linefit:
    def __init__(self):
        self.e0 = None
        self.e1 = None
        self.cbin_idx = None
        self.cbin = None
        self.res_arr = None
        self.smax = None
        self.step_t = None
        self.t_ini = None


    def get_linefit_static_data(self):
        if self.cbin == self.smax:
            self.nbin = self.cbin
            self.bin_idx = self.cbin_idx
        else:
            e_list = [self.e0, self.e1]
            eout, cBIdx, cbin = cut_resolution_for_linefit(
                e_list, self.cbin_idx, self.res_arr, self.smax
            )
            self.e0, self.e1 = [eout[0, :, :, :], eout[1, :, :, :]]
            self.bin_idx = cBIdx
            self.nbin = cbin

    def get_fsc_wght(self, e0, ert, bin_idx, nbin):
        cx, cy, cz = e0.shape
        bin_stats = core.fsc.anytwomaps_fsc_covariance(e0, ert, bin_idx, nbin)
        fsc, _ = bin_stats[0], bin_stats[1]
        fsc = np.array(fsc, dtype=np.float64, copy=False)
        w_grid = fcodes_fast.read_into_grid(bin_idx, fsc, nbin, cx, cy, cz)
        w_grid = np.array(w_grid, dtype=np.float64, copy=False)
        return w_grid

    def func_t(self, i):
        nx, ny, nz = self.e0.shape
        t = self.step_t * i
        st, _, _, _ = fcodes_fast.get_st(nx, ny, nz, t)
        e1_t = self.e1 * st
        w_grid = self.get_fsc_wght(self.e0, e1_t, self.bin_idx, self.nbin)
        fval = np.sum(w_grid * self.e0 * np.conjugate(e1_t))
        return -fval.real

    def scalar_opt_trans(self):
        from scipy.optimize import minimize_scalar

        f = self.func_t
        res = minimize_scalar(f, method="brent")
        return res.x


class EmFit:
    def __init__(self, mapobj, interp="linear", dfs=None):
        self.mapobj = mapobj
        self.cut_dim = mapobj.cdim
        self.ful_dim = mapobj.map_dim
        self.cell = mapobj.map_unit_cell
        self.pixsize = mapobj.pixsize
        self.origin = mapobj.map_origin
        #self.bin_idx = mapobj.bin_idx
        #self.nbin = mapobj.nbin
        self.interp = interp
        self.dfs = dfs
        self.w_grid = None
        self.fsc = None
        self.sv = None
        self.t = None
        self.st = None
        self.step = None
        self.q = None
        self.q_accum = None
        self.q_final_list = []
        self.rotmat = None
        self.t_accum = None
        self.ert = None
        self.frt = None
        self.e0 = None
        self.e1 = None
        self.w2_grid = None
        self.fsc_lst = []

    def calc_fsc(self):
        cx, cy, cz = self.e0.shape
        self.st, s1, s2, s3 = fcodes_fast.get_st(cx, cy, cz, self.t)
        self.sv = np.array([s1, s2, s3])
        self.ert = self.e1 * self.st
        fsc = core.fsc.anytwomaps_fsc_covariance(
            self.e0, self.ert, self.mapobj.cbin_idx, self.mapobj.cbin
        )[0]
        return fsc

    def get_wght(self):
        cx, cy, cz = self.e0.shape
        w_grid = fcodes_fast.read_into_grid(
            self.mapobj.cbin_idx,
            self.fsc, #/ (1 - self.fsc ** 2),
            self.mapobj.cbin,
            cx,
            cy,
            cz,
        )
        fsc_sqd = self.fsc ** 2
        fsc_combi = fsc_sqd #/ (1 - fsc_sqd)
        w2_grid = fcodes_fast.read_into_grid(
            self.mapobj.cbin_idx, fsc_combi, self.mapobj.cbin, cx, cy, cz
        )
        return w_grid, w2_grid

    def functional(self):
        fval = np.sum(self.w_grid * self.e0 * np.conjugate(self.ert))
        return fval.real

    def minimizer(self, ncycles, t_init, rotmat, ifit, smax_lf, tol=1.e-2, fobj=None, q_init=None):
        fsc_lst = []
        fval_list = []
        q_list = []
        t_list = []
        self.e0 = self.mapobj.ceo_lst[0]  # Static map e-data for fit
        self.e1 = self.mapobj.ceo_lst[1]
        print("pixel size: ", self.pixsize)
        print("Cycle#   ", "Fval  ", "  Trans(A)  ", "  avg(FSC)  ")
        for i in range(ncycles):
            start = timer()
            if i == 0:
                self.t = np.asarray(t_init, dtype="float")
                t_accum = self.t
                translation_vec = trans2angstroms(t_accum, self.pixsize, self.mapobj.map_dim)
            # check FSC and return parameters accordingly
            self.fsc = self.calc_fsc()
            self.w_grid, self.w2_grid = self.get_wght()
            fval = self.functional()
            fval_list.append(fval)
            t_list.append(t_accum)
            if i == 0:
                fval_previous = fval
                fsc = self.fsc
                fsc_lst.append(fsc)
            if i > 0 and fval_previous < fval or i == ncycles - 1:
                fsc = self.fsc
            if i > 0 and i == ncycles - 1:
                self.t_accum = t_accum
                self.fsc_lst = fsc_lst
                translation_vec = trans2angstroms(self.t_accum, self.pixsize, self.mapobj.map_dim)
                #t_accum_angstrom = self.t_accum * self.pixsize[0] * self.mapobj.map_dim[0]
                #translation_vec = np.sqrt(np.dot(t_accum_angstrom, t_accum_angstrom))
                break
            if i > 0 and abs(fval_list[-1] - fval_list[-2]) < tol:
                self.t_accum = t_accum
                break
            print(
                "{:5d} {:8.4f} {:6.4f} {:6.4f}".format(
                    i, fval, translation_vec, np.average(self.fsc)
                )
            )
            self.step = derivatives_translation(self.e0, self.ert, self.w_grid, self.w2_grid, self.sv)
            self.e1 = self.ert
            lft = linefit()
            lft.cbin_idx = self.mapobj.cbin_idx
            lft.cbin = self.mapobj.cbin
            lft.res_arr = self.mapobj.res_arr
            lft.smax = smax_lf
            lft.e0 = self.e0
            lft.e1 = self.e1
            lft.get_linefit_static_data()
            lft.step_t = self.step
            alpha = lft.scalar_opt_trans()
            # translation
            self.t = self.step[:3] * alpha
            t_accum = t_accum + self.t
            t_accum_angstrom = t_accum * self.pixsize[0] * self.mapobj.map_dim[0]
            translation_vec = np.sqrt(np.dot(t_accum_angstrom, t_accum_angstrom))
            fval_previous = fval
            end = timer()
            if timeit:
                print("time for one cycle:", end - start)


def fsc_between_static_and_transfomed_map(
        staticmap, movingmap, bin_idx, t, nbin):
    nx, ny, nz = staticmap.shape
    st, _, _, _ = fcodes_fast.get_st(nx, ny, nz, t)
    f1f2_fsc = core.fsc.anytwomaps_fsc_covariance(staticmap, 
        movingmap * st, bin_idx, nbin)[0]
    return f1f2_fsc


def run_fit(
    emmap1,
    #rotmat,
    t,
    ncycles,
    ifit,
    fobj=None,
    interp=None,
    fitres=None,
    dfs_full=None,
):
    ifit = 1
    if fitres is not None:
        if fitres <= emmap1.res_arr[-1]:
            fitbin = len(emmap1.res_arr) - 1
        else:
            dist = np.sqrt((emmap1.res_arr - fitres) ** 2)
            ibin = np.argmin(dist)
            if ibin % 2 != 0:
                ibin = ibin - 1
            fitbin = min([len(dist), ibin])
    if fitres is None:
        fitbin = len(emmap1.res_arr) - 1
    fsc_lst = []

    q_final = np.array([1.0, 0.0, 0.0, 0.0], dtype=np.float64)
    rotmat = np.identity(3)
    import emda.ext.utils as utils
    for i in range(5):
        print("Resolution cycle #: ", i)
        if i == 0:
            f1f2_fsc = fsc_between_static_and_transfomed_map(
                staticmap=emmap1.fo_lst[0],
                movingmap=emmap1.fo_lst[ifit],
                bin_idx=emmap1.bin_idx,
                t=t,
                nbin=emmap1.nbin,
            )
            #print('starting FSC: ')
            #print(f1f2_fsc)
            fsc_lst.append(f1f2_fsc)
            ibin = utils.determine_ibin(f1f2_fsc)
            if fitbin < ibin:
                ibin = fitbin
            print("Fitting starts at ", emmap1.res_arr[ibin], " (A)")
            ibin_old = ibin
            if np.average(f1f2_fsc) > 0.999:
                rotmat = rotmat
                t = t
                fsc_lst.append(f1f2_fsc)
                print("\n***FSC between static and moving maps***\n")
                print("bin#     resolution(A)      start-FSC     end-FSC\n")
                for j in range(len(emmap1.res_arr)):
                    print(
                        "{:5d} {:6.2f} {:8.4f} {:8.4f}".format(
                            j, emmap1.res_arr[j], fsc_lst[1][j], fsc_lst[0][j]
                        )
                    )
                break
        else:
            # Apply translation to calculate current fsc
            f1f2_fsc = fsc_between_static_and_transfomed_map(
                emmap1.fo_lst[0],
                emmap1.fo_lst[ifit],
                emmap1.bin_idx,
                t,
                emmap1.nbin,
            )
            ibin = utils.determine_ibin(f1f2_fsc)
            if fitbin < ibin:
                ibin = fitbin
            if ibin_old == ibin:
                fsc_lst.append(f1f2_fsc)
                print("\n***FSC between static and moving maps***\n")
                print("bin#     resolution(A)      start-FSC     end-FSC\n")
                for j in range(len(emmap1.res_arr)):
                    print(
                        "{:5d} {:6.2f} {:8.4f} {:8.4f}".format(
                            j, emmap1.res_arr[j], fsc_lst[0][j], fsc_lst[1][j]
                        )
                    )
                break
            else:
                ibin_old = ibin
        if ibin == 0:
            print("ibin = 0")
            raise SystemExit("Cannot proceed! Stopping now...")

        e_list = [emmap1.eo_lst[0], emmap1.eo_lst[ifit]]
        eout, cBIdx, cbin = cut_resolution_for_linefit(
            e_list, emmap1.bin_idx, emmap1.res_arr, ibin
        )
        emmap1.ceo_lst = [eout[0, :, :, :], eout[1, :, :, :]]
        emmap1.cdim = eout[1, :, :, :].shape
        emmap1.cbin_idx = cBIdx
        emmap1.cbin = cbin
        rfit = EmFit(emmap1)
        slf = min([ibin, 50])
        rfit.minimizer(ncycles, t, rotmat, ifit, smax_lf=slf, fobj=fobj)
        t = rfit.t_accum
    return t, q_final


def overlay(
    maplist,
    ncycles=50,
    t_init=[0.0, 0.0, 0.0],
    symorder=3,
    rotaxis=[0,0,1],
    interp='linear',
    fitres=4,
):
    emmap1 = EmmapOverlay(maplist)
    emmap1.load_maps()
    emmap1.calc_fsc_from_maps()
    t = [itm / emmap1.pixsize[i] for i, itm in enumerate(t_init)]
    trans_list = []
    
    # Case 1: 
    # There are two map and 
    # we want to optimize translation of the 2nd map to 1st map
    """ for ifit in range(1, len(emmap1.eo_lst)):
        t, _ = run_fit(
            emmap1=emmap1,
            #rotmat=rotmat_init,
            t=t,
            ncycles=ncycles,
            ifit=ifit,
            interp=interp,
            fitres=fitres,
        )
        trans_list.append(t) """

    # Case 2: 
    # We have one map and known rotation (axis and rotation angle) 
    # and we want to optimize the translation between the original map
    # and its symmetry copies
    axis = np.asarray(rotaxis)
    axis = axis / math.sqrt(np.dot(axis, axis))
    q_init_list = []
    for n in range(1, symorder):
        theta = float(n * 360. / symorder)
        q = core.quaternions.get_quaternion([[axis[2], axis[1], axis[0]], [theta]])
        print(theta, q)
        q_init_list.append(q)
    f0 = emmap1.fo_lst[0]
    for ifit, q in enumerate(q_init_list):
        rotmat = core.quaternions.get_RM(q)
        # apply rotation
        frs = fcodes_fast.trilinear2(f0,
                                     bin_idx=emmap1.bin_idx,
                                     rm=rotmat,
                                     nbin=emmap1.nbin,
                                     mode=0,
                                     ncopies=1,
                                     nz=f0.shape[2],
                                     ny=f0.shape[1],
                                     nx=f0.shape[0])[:,:,:,0]
        emmap1.fo_lst = [f0, frs]
        emmap1.eo_lst = [f0, frs]
        t, _ = run_fit(
            emmap1=emmap1,
            #rotmat=rotmat_init,
            t=t,
            ncycles=ncycles,
            ifit=ifit,
            interp=interp,
            fitres=fitres,
        )
        trans_list.append(t)
    # output maps
    output_rotated_maps(emmap1, trans_list)
    return emmap1, trans_list


def output_rotated_maps(emmap1, tlist):
    cell = emmap1.map_unit_cell
    bin_idx = emmap1.bin_idx
    nbin = emmap1.nbin
    f0 = emmap1.fo_lst[0]
    f1 = emmap1.fo_lst[1]
    nz, ny, nx = f0.shape
    data2write = np.real(ifftshift(ifftn(ifftshift(f0))))
    core.iotools.write_mrc(data2write, "static_map.mrc", cell)
    f1f2_fsc_unaligned = core.fsc.anytwomaps_fsc_covariance(
            f0, f1, bin_idx, nbin)[0]
    # average t
    print('tlist: ', tlist)
    t = np.mean(np.asarray(tlist), axis=0)
    print('average translation vector (A): ')
    print(t * np.asarray(emmap1.pixsize, 'float') * np.asarray(emmap1.map_dim, 'int'))

    st, _, _, _ = fcodes_fast.get_st(nz, ny, nx, t)
    frt = f1 * st
    data2write = np.real(ifftshift(ifftn(ifftshift(frt))))
    core.iotools.write_mrc(data2write,"fitted_map.mrc",cell)
    f1f2_fsc = core.fsc.anytwomaps_fsc_covariance(
            f0, frt, bin_idx, nbin)[0]
    core.plotter.plot_nlines(
            emmap1.res_arr,
            [f1f2_fsc_unaligned[:emmap1.nbin], f1f2_fsc[:emmap1.nbin]],
            "fsc.eps",["FSC before", "FSC after"])
    


if __name__ == "__main__":
    maplist = [
        #"/Users/ranganaw/MRC/REFMAC/haemoglobin/EMD-3651/emda_test/map_transform/emd_3651.map",
        #"/Users/ranganaw/MRC/REFMAC/haemoglobin/EMD-3651/emda_test/map_transform/transformed.mrc"
        #"/Users/ranganaw/MRC/REFMAC/EMD-6952/map/emd_6952.map"
        "/Users/ranganaw/MRC/REFMAC/Keitaro/EMD-21951/emd_21951.map"
        ] 
    
    emmap1, transl_lst = overlay(maplist,
                                 ncycles=50,
                                 t_init=[0.0, 0.0, 0.0],
                                 symorder=2,
                                 rotaxis=[1,0,0],
                                 fitres=6,
                                 )
 