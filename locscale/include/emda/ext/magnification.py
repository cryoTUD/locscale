"""
Author: "Rangana Warshamanage, Garib N. Murshudov"
MRC Laboratory of Molecular Biology
    
This software is released under the
Mozilla Public License, version 2.0; see LICENSE.
"""

from __future__ import absolute_import, division, print_function, unicode_literals
import numpy as np
from emda.core import iotools, restools
import fcodes_fast as fcodes


class LineFit:
    def __init__(self):
        self.e0_lf = None
        self.step = None
        self.e1_lf = None
        self.cbin_idx_lf = None
        self.cbin_lf = None
        self.w_grid = None
        self.res_arr = None
        self.mode = "model"

    def get_linefit_static_data(self, e0, cbin_idx, res_arr, smax):
        (
            self.e0_lf,
            self.cbin_idx_lf,
            self.cbin_lf,
        ) = restools.cut_resolution_for_linefit(e0, cbin_idx, res_arr, smax)
        self.res_arr = res_arr[: self.cbin_lf]

    def f(self, kini):
        import emda.ext.mapfit.utils as utils

        nx, ny, nz = self.e0_lf.shape
        ncopies = 1
        k = self.k + kini[0] * self.step
        eck = fcodes.tricubic_zoom(
            float(k), self.e1_lf, 0, ncopies, nx, ny, nz)
        eckt = eck[:, :, :, 0]
        w_grid = utils.get_fsc_wght(
            self.e0_lf, eckt, self.cbin_idx_lf, self.cbin_lf)
        fval = np.sum(w_grid * self.e0_lf * np.conjugate(eckt))
        return -fval.real

    def fn2(self, kini):
        nx, ny, nz = self.e0_lf.shape
        ncopies = 1
        #k = self.k + kini[0] * self.step
        k = self.k + kini * self.step
        fcks = fcodes.tricubic_zoom(
            float(k), self.e1_lf, 0, ncopies, nx, ny, nz)
        fckt = fcks[:, :, :, 0]
        _, _, _, fval = fcodes.scalesigmafval_full(
            self.e0_lf,
            fckt,
            self.cbin_idx_lf,
            self.res_arr,
            0,
            self.cbin_lf,
            nx,
            ny,
            nz,
        )
        return fval.real

    def calc_fval_for_different_kvalues_at_this_step(self, k=1.0, e1=None):
        from scipy import optimize

        # start = timer()
        nx, ny, nz = e1.shape
        cx = self.e0_lf.shape[0] // 2
        cy = cz = cx
        dx = int((nx - 2 * cx) / 2)
        dy = int((ny - 2 * cy) / 2)
        dz = int((nz - 2 * cz) / 2)
        self.e1_lf = e1[dx: dx + 2 * cx, dy: dy + 2 * cy, dz: dz + 2 * cz]
        assert self.e1_lf.shape == self.e0_lf.shape
        init_guess = [k]
        self.k = k
        if self.mode == "map":  # for map-map
            minimum = optimize.minimize(self.f, init_guess, method="Powell")
        if self.mode == "model":  # for map-model
            #minimum = optimize.minimize(self.fn2, init_guess, method="Powell")
            minimum = optimize.minimize_scalar(self.fn2, method="brent")
        # end = timer()
        # print(' time for line search: ', end-start)
        return minimum.x


def create_xyz_grid(uc, nxyz):
    x = np.fft.fftfreq(nxyz[0])  # * uc[0]
    y = np.fft.fftfreq(nxyz[1])  # * uc[1]
    z = np.fft.fftfreq(nxyz[2])  # * uc[2]
    xv, yv, zv = np.meshgrid(x, y, z)
    xyz = [yv, xv, zv]
    for i in range(3):
        xyz[i] = np.fft.ifftshift(xyz[i])
    return xyz


def get_xyz_sum(xyz):
    xyz_sum = np.zeros(shape=(6), dtype="float")
    n = -1
    for i in range(3):
        for j in range(3):
            if i == 0:
                sumxyz = np.sum(xyz[i] * xyz[j])
            elif i > 0 and j >= i:
                sumxyz = np.sum(xyz[i] * xyz[j])
            else:
                continue
            n = n + 1
            xyz_sum[n] = sumxyz
    return xyz_sum


def dFks(mapin, uc):
    xyz = create_xyz_grid(uc, mapin.shape)
    vol = uc[0] * uc[1] * uc[2]

    # Calculating dFC(ks)/dk using FFT
    dfk = np.zeros(mapin.shape, np.complex64)
    xyz_sum = 0.0
    for i in range(3):
        xyz_sum = xyz_sum + np.sum(xyz[i])
    # dfk = np.fft.fftshift((-1/vol) * 2j * np.pi * np.fft.fftn(mapin * xyz_sum))
    dfk = np.fft.fftshift(-2j * np.pi * np.fft.fftn(mapin * xyz_sum))

    # second derivative
    xyz_sum = 0.0
    tp2 = (2.0 * np.pi) ** 2
    """ for i in range(3):
        for j in range(3):
            xyz_sum = xyz_sum + np.sum(xyz[i] * xyz[j]) """
    xyz_sum = np.sum(get_xyz_sum(xyz))
    # ddfk = -(tp2/vol) * np.fft.fftshift(np.fft.fftn(mapin * xyz_sum))
    ddfk = -(tp2) * np.fft.fftshift(np.fft.fftn(mapin * xyz_sum))
    return dfk, ddfk


def dFts(Fc, sv):
    nx, ny, nz = Fc.shape
    tp2 = (2.0 * np.pi) ** 2
    dt_arr = np.zeros(shape=(nx, ny, nz, 3), dtype="complex")
    ddt_arr = np.zeros(shape=(nx, ny, nz, 3, 3), dtype="complex")
    for i in range(3):
        # 1st derivative
        dfs = 2.0 * 1j * np.pi * sv[i] * Fc
        dt_arr[:, :, :, i] = dfs
        # 2nd derivative
        for j in range(3):
            if i == 0:
                ddfs = -tp2 * sv[i] * sv[j] * Fc
            elif i > 0 and j >= i:
                ddfs = -tp2 * sv[i] * sv[j] * Fc
            else:
                ddfs = ddt_arr[:, :, :, j, i]
            ddt_arr[:, :, :, i, j] = ddfs
    return dt_arr, ddt_arr


def get_ll(emmap1, fmaplist, bin_idx, res_arr, nbin, k, smax):
    fref = fmaplist[0]
    ftar = fmaplist[1]
    nx, ny, nz = fref.shape  # model
    ncopies = 1
    fcks = fcodes.tricubic_zoom(k, fref, 0, ncopies, nx, ny, nz)
    # optimize translation
    t = np.array([0.0, 0.0, 0.0], 'float')
    f1 = fcks[:, :, :, 0]
    t = optimise_translation(f1, ftar, t, bin_idx, nbin, res_arr, smax)
    t = np.asarray(t)
    st, s1, s2, s3 = fcodes.get_st(nx, ny, nz, t)
    sv = np.array([s1, s2, s3])
    #
    fckt = fcks[:, :, :, 0] * st
    if len(fmaplist) == 2:
        scale_d, sigma, totalvar, fval = fcodes.scalesigmafval_full(
            ftar, fckt, bin_idx, res_arr, 0, nbin, nx, ny, nz
        )
        return ftar, fckt, scale_d, sigma, totalvar, fval, sv, t


def derivatives_mapmodel(fo, fc, bin_idx, sv, D, totalvar, uc):
    # 1. Calculate dFc/dk and dFc/dT
    mapin = np.fft.ifftn(np.fft.ifftshift(fc))
    dk, ddk = dFks(mapin, uc)
    # dt, ddt = dFts(fc,sv)

    nbin = len(totalvar)
    nx, ny, nz = fc.shape
    dll, ddll = fcodes.ll_derivatives(
        fo, fc, bin_idx, D, totalvar, 1, nbin, nx, ny, nz)

    # 1st derivatives
    df_val = np.zeros(shape=(4), dtype="float")
    df_val[0] = np.sum(np.real(dll * np.conjugate(dk)))
    """ for i in range(3):
        df_val[i+1] = np.sum(np.real(dll * np.conjugate(2j * np.pi * sv[i] * fc))) """

    # 2nd derivatives
    ddf_val = np.zeros(shape=(4, 4), dtype="float")
    ddf_val[0, 0] = np.sum(np.real(ddll * np.conjugate(ddk)))
    """ tp2 = (2.0 * np.pi)**2
    for i in range(3):
        for j in range(3):
            ddf_val[i+1,j+1] = -tp2 * np.sum(np.real(ddll * \
                np.conjugate(fc * sv[i] * sv[j]))) """
    ddf_val_inv = np.linalg.pinv(ddf_val)
    step = ddf_val_inv.dot(-df_val)
    return step


def calc_fsc_t(fstatic, frotated, t, bin_idx, nbin):
    from emda import core

    cx, cy, cz = fstatic.shape
    st, s1, s2, s3 = fcodes.get_st(cx, cy, cz, t)
    sv = np.array([s1, s2, s3])
    frt = frotated * st
    fsc = core.fsc.anytwomaps_fsc_covariance(fstatic, frt, bin_idx, nbin)[
        0
    ]
    return frt, fsc, sv


def get_wght(f0, fsc, bin_idx, nbin):
    cx, cy, cz = f0.shape
    w_grid = fcodes.read_into_grid(
        bin_idx,
        fsc,  # fsc / (1 - fsc ** 2),
        nbin,
        cx,
        cy,
        cz,
    )
    fsc_sqd = fsc ** 2
    fsc_combi = fsc_sqd,  # fsc_sqd / (1 - fsc_sqd)
    w2_grid = fcodes.read_into_grid(
        bin_idx, fsc_combi, nbin, cx, cy, cz
    )
    return w_grid, w2_grid


def optimise_translation(f0, f1, t, bin_idx, nbin, res_arr, smax_lf, ncy=2):
    from emda.ext.sym.refine_symaxis import derivatives_translation, linefit2

    #print("translation optimisation..")
    for i in range(ncy):
        frt, fsc, sv = calc_fsc_t(f0, f1, t, bin_idx, nbin)
        w_grid, w2_grid = get_wght(f0, fsc, bin_idx, nbin)
        #fval = np.sum(w_grid * f0 * np.conjugate(frt))
        #print(i, fval.real, t)
        step = derivatives_translation(f0, frt, w_grid, w2_grid, sv)
        lft = linefit2()
        lft.get_linefit_static_data(
            [f0, frt], bin_idx, res_arr, smax_lf
        )
        lft.step = step
        alpha = lft.scalar_opt_trans()
        t = t + step * alpha
    return t


def apply_transformation(folist, rmlist):
    assert len(rmlist) == len(folist)
    #print(len(folist))
    #print(folist[0].shape)
    nx, ny, nz = folist[0].shape
    frs = fcodes.trilinearn(np.stack((folist[0], folist[1]), axis = 0),
                            np.stack((rmlist[0], rmlist[1]), axis = 0),
                            0,
                            len(folist),
                            nx, ny, nz)
    #print(frs.shape)
    return frs


def apply_translation(fo, t):
    nx, ny, nz = fo.shape
    st, _, _, _ = fcodes.get_st(nx, ny, nz, t)    
    return fo * st


def cut_resolution_for_linefit(farr, bin_idx, res_arr, smax):
    # Making data for map fitting
    ncopies, nx, ny, nz = farr.shape
    cbin = cx = smax
    dx = int((nx - 2 * cx) / 2)
    dy = int((ny - 2 * cx) / 2)
    dz = int((nz - 2 * cx) / 2)
    cBIdx = bin_idx[dx : dx + 2 * cx, dy : dy + 2 * cx, dz : dz + 2 * cx]
    fout = fcodes.cutmap_arr(
        farr, bin_idx, cbin, 0, len(res_arr), nx, ny, nz, ncopies
    )[:, dx : dx + 2 * cx, dy : dy + 2 * cx, dz : dz + 2 * cx]
    return fout, cBIdx, cbin


def minimizer_mapmodel(emmap1, resol=4.0, ncycles=10, rmlist=None, tlist=None):
    print('\nOptimising map magnification...')
    tol = 1e-2 # tolerence for refinement convergence
    # apply transformation
    frs = apply_transformation(folist=emmap1.fo_lst, rmlist=rmlist)
    # output reference map
    refmap = np.fft.fftshift(
        (np.fft.ifftn(np.fft.ifftshift(frs[0,:,:,:]))).real)
    iotools.write_mrc(refmap, 'reference.mrc', emmap1.map_unit_cell)

    # find the bin for given resolution for fit
    dist = np.sqrt((emmap1.res_arr[:emmap1.nbin] - resol) ** 2)
    smax = np.argmin(dist)
    # cut data
    fout, cBIdx, cbin = cut_resolution_for_linefit(
        frs, emmap1.bin_idx, emmap1.res_arr, smax)

    # parameters for magnification refinement
    bin_idx = cBIdx
    res_arr = emmap1.res_arr[:cbin]
    nbin = cbin
    uc = emmap1.map_unit_cell

    for ifit in range(1, fout.shape[0]):
        # output moving map before magref refinement
        movingmap = apply_translation(frs[ifit,:,:,:], tlist[ifit])
        refmap = np.fft.fftshift(
                    (np.fft.ifftn(np.fft.ifftshift(movingmap))).real)
        iotools.write_mrc(refmap, 'startmap_'+str(ifit)+'.mrc', uc)
        # use cutdata in mag. refinement
        maplist = [fout[0,:,:,:], apply_translation(fout[ifit,:,:,:], tlist[ifit])]
        # initial parameters
        t = [0.0, 0.0, 0.0]
        k = k_previous = 1.0
        for i in range(ncycles):
            fo, fckt, scale_d, sigma, totalvar, fval, sv, t = get_ll(
                emmap1, maplist, bin_idx, res_arr, nbin, k, smax)
            if i == 0:
                fval_previous = fval
                print()
                print("ifit    cycle#    func val.   magnification")
                print(
                    "{:5d} {:5d} {:8.4f} {:6.4f}".format(
                        ifit, i, fval, k
                    )
                )
            if i > 0 and abs(fval_previous - fval) > tol or i == ncycles - 1:
                print(
                    "{:5d} {:5d} {:8.4f} {:6.4f}".format(
                        ifit, i, fval, k
                    )
                )
                k_previous = k
            if i > 0 and abs(fval - fval_previous) <= tol or i == ncycles - 1:
                print(
                    "{:5d} {:5d} {:8.4f} {:6.4f}".format(
                        ifit, i, fval, k
                    )
                )
                # output map after mag. refinement
                magerror = abs(k_previous - 1.0) * 100.0
                print("magnification error (%): ", "{:.3f}".format(magerror))
                _, nx, ny, nz = frs.shape
                fck = fcodes.tricubic_zoom(
                    float(1 / k_previous), frs[ifit,:,:,:], 0, 1, nx, ny, nz)[:, :, :, 0]
                t = optimise_translation(frs[0,:,:,:], 
                                         fck, 
                                         t, 
                                         emmap1.bin_idx, 
                                         emmap1.nbin, 
                                         emmap1.res_arr, 
                                         smax, 4)
                fckt = fck * fcodes.get_st(nx, ny, nz, t)[0]
                zoomedmap = np.fft.fftshift(
                    (np.fft.ifftn(np.fft.ifftshift(fckt))).real)
                mapname = "emda_magcorretedmap_" + str(ifit) + ".mrc"
                iotools.write_mrc(zoomedmap, mapname, uc)
                print("Magnification corrected map %s was written." %(mapname) )
                break
            step = derivatives_mapmodel(
                fo, fckt, bin_idx, sv, scale_d, sigma, uc)
            if i == 0:
                linefit = LineFit()
                linefit.get_linefit_static_data(
                    e0=fo, cbin_idx=bin_idx, res_arr=res_arr, smax=smax
                )
            linefit.step = step[0]
            alpha = linefit.calc_fval_for_different_kvalues_at_this_step(
                k=k, e1=maplist[0])
            k = k + alpha * step[0]
            fval_previous = fval


def prepare_data(maplist, masklist=None):
    from emda.ext.overlay import EmmapOverlay
    try:
        emmap1 = EmmapOverlay(maplist, masklist)
    except:
        emmap1 = EmmapOverlay(maplist)
    emmap1.com = False
    emmap1.load_maps()
    emmap1.calc_fsc_from_maps()
    return emmap1


def optimize_superposition(emmap1, ncycles=50):
    from emda.ext.overlay import run_fit
    from emda.core.quaternions import get_RM

    print("\nOptimising overlay.....\n")
    rotmat = np.identity(3)
    t = [0.0, 0.0, 0.0]
    rotmat_list = []
    trans_list = []
    fobj = open("EMDA_overlay.txt", "w")
    rotmat_list.append(rotmat)
    trans_list.append(t)
    for ifit in range(1, len(emmap1.eo_lst)):
        t, q_final = run_fit(
                emmap1,
                rotmat,
                t,
                ncycles,
                ifit,
            )
        rotmat_list.append(get_RM(q_final))
        trans_list.append(t)
    print("Optimising overlay.....Done")
    return rotmat_list, trans_list

def main(maplist, masklist=None, resol=None):
    emmap1 = prepare_data(maplist, masklist)
    rm_list, t_list = optimize_superposition(emmap1)
    minimizer_mapmodel(emmap1=emmap1, resol=resol, rmlist=rm_list, tlist=t_list)


if __name__ == "__main__":
    maplist = ["/Users/ranganaw/MRC/REFMAC/haemoglobin/EMD-3651/other/emd_3651_half_map_1.map",
                "/Users/ranganaw/MRC/REFMAC/haemoglobin/EMD-3651/emda_test/magnification/data_for_emda-paper/emd_3651_half1_k095.mrc"]
    #maplist = [
    #    "/Users/ranganaw/MRC/REFMAC/beta_gal/magref_for_paper/xtallography_reference.mrc",
    #    "/Users/ranganaw/MRC/REFMAC/beta_gal/magref_for_paper/emd_10574_fullmap_resampled2staticmap_pca.map",
    #           ]
    main(maplist)
