# lock rotation
from __future__ import absolute_import, division, print_function, unicode_literals
from timeit import default_timer as timer
import numpy as np
import sys, math
import fcodes_fast
import emda.emda_methods as em
from numpy.fft import fftn, ifftn, fftshift, ifftshift
from emda.ext.mapfit import map_class
from emda import core
from emda.core import quaternions
from emda.ext.mapfit.utils import (
    get_FRS,
    create_xyz_grid,
    get_xyz_sum,
    set_dim_even,
)
from emda.core.quaternions import derivatives_wrt_q
from emda.ext.mapfit.derivatives import get_dqda, new_dFs2
from scipy.optimize import minimize_scalar
from matplotlib import pyplot as plt
from emda.ext.mapfit import frequency_marching
from emda.ext.sym import run_proshade as proshade
from more_itertools import sort_together
#from rebox_map import reboxmap

mapout = True

class EmmapOverlay:
    def __init__(self, imap, imask=None):
        self.imap = imap
        self.imask = imask
        self.map_unit_cell = None
        self.map_origin = None
        self.map_dim = None
        self.pixsize = None
        self.arr_lst = []
        self.arr_original = None
        self.ceo_lst = None
        self.cfo_lst = None
        self.cbin_idx = None
        self.cdim = None
        self.cbin = None
        self.com = True
        self.com1 = None
        self.box_centr = None
        self.fhf_lst = []
        self.nbin = None
        self.res_arr = None
        self.bin_idx = None
        self.fo_lst = None
        self.eo_lst = None
        self.totalvar_lst = None
        self.q_init_list = []

    def get_maps(self):
        uc, arr, origin = em.get_data(self.imap)
        arr = set_dim_even(arr)
        self.arr_original = arr
        com = em.center_of_mass_density(arr)
        print("com:", com)
        nx, ny, nz = arr.shape
        box_centr = (nx // 2, ny // 2, nz // 2)
        self.box_centr = box_centr
        self.com1 = com
        if self.com:
            arr = em.shift_density(arr, np.subtract(box_centr, com))
            print("com after centering:", em.center_of_mass_density(arr))
        self.arr_lst.append(arr)
        self.fhf_lst.append(fftshift(fftn(fftshift(arr))))
        self.pixsize = uc[0] / arr.shape[0]
        self.map_origin = origin
        self.map_unit_cell = uc
        self.map_dim = arr.shape

    def calc_fsc_from_maps(self):
        # function for only two maps fitting
        nmaps = len(self.fhf_lst)
        fFo_lst = []
        fEo_lst = []
        fBTV_lst = []
        self.nbin, self.res_arr, self.bin_idx = core.restools.get_resolution_array(
            self.map_unit_cell, self.fhf_lst[0]
        )
        for i in range(nmaps):
            _, _, _, totalvar, fo, eo = core.fsc.halfmaps_fsc_variance(
                self.fhf_lst[i], self.fhf_lst[i], self.bin_idx, self.nbin
            )
            fFo_lst.append(fo)
            fEo_lst.append(eo)
            fBTV_lst.append(totalvar)
        self.fo_lst = fFo_lst
        self.eo_lst = fEo_lst
        self.totalvar_lst = fBTV_lst


def cut_resolution_for_linefit(f_list, bin_idx, res_arr, smax):
    # Making data for map fitting
    f_arr = np.asarray(f_list, dtype="complex")
    nx, ny, nz = f_list[0].shape
    cbin = cx = smax
    dx = int((nx - 2 * cx) / 2)
    dy = int((ny - 2 * cx) / 2)
    dz = int((nz - 2 * cx) / 2)
    cBIdx = bin_idx[dx: dx + 2 * cx, dy: dy + 2 * cx, dz: dz + 2 * cx]
    fout = fcodes_fast.cutmap_arr(
        f_arr, bin_idx, cbin, 0, len(res_arr), nx, ny, nz, len(f_list)
    )[:, dx: dx + 2 * cx, dy: dy + 2 * cx, dz: dz + 2 * cx]
    return fout, cBIdx, cbin


class linefit:
    def __init__(self):
        self.e0 = None
        self.e1 = None
        self.bin_idx = None
        self.nbin = None
        self.axis = None
        self.angle_list = []
        self.step_axis = None
        self.step_t = None
        self.t_ini = None
        self.alpha_t = None
        self.axis_ini = None

    def get_linefit_static_data(self, e_list, bin_idx, res_arr, smax):
        if len(e_list) == 2:
            eout, self.bin_idx, self.nbin = cut_resolution_for_linefit(
                e_list, bin_idx, res_arr, smax
            )
        else:
            print("len(e_list: ", len(e_list))
            raise SystemExit()
        self.e0 = eout[0, :, :, :]
        self.e1 = eout[1, :, :, :]

    def get_quaternion(self, axis, angle_i):
        rv = axis
        s = math.sin(angle_i / 2.0)
        q1, q2, q3 = rv[0] * s, rv[1] * s, rv[2] * s
        q0 = math.cos(angle_i / 2.0)
        q = np.array([q0, q1, q2, q3], dtype=np.float64)
        return q

    def get_fsc_wght(self, e0, ert, bin_idx, nbin):
        cx, cy, cz = e0.shape
        bin_stats = core.fsc.anytwomaps_fsc_covariance(e0, ert, bin_idx, nbin)
        fsc, _ = bin_stats[0], bin_stats[1]
        fsc = np.array(fsc, dtype=np.float64, copy=False)
        w_grid = fcodes_fast.read_into_grid(bin_idx, fsc, nbin, cx, cy, cz)
        w_grid = np.array(w_grid, dtype=np.float64, copy=False)
        return w_grid

    def func(self, i):
        step_ax = self.step_axis * i
        self.axis = getaxis(self.axis_ini, step_ax)
        fval = 0.0
        q_i = self.get_quaternion(self.axis, self.angle_list[0])
        rotmat = core.quaternions.get_RM(q_i)
        ers = get_FRS(rotmat, self.e0, interp="linear")
        w_grid = self.get_fsc_wght(
            self.e0, ers[:, :, :, 0], self.bin_idx, self.nbin)
        fval = np.real(np.sum(w_grid * self.e0 *
                              np.conjugate(ers[:, :, :, 0])))
        return -fval

    def scalar_opt(self, t=None):
        f = self.func
        res = minimize_scalar(f, method="brent")
        return res.x


def derivatives(e0, e1, w_grid, w2_grid, q, sv, xyz, xyz_sum, vol):


    nx, ny, nz = e0.shape
    sv_np = np.zeros((nx, ny, nz, 3), dtype=np.float64)
    for i in range(3):
        sv_np[:, :, :, i] = sv[i]
    dRdq = derivatives_wrt_q(q)
    dqda = get_dqda(q)
    dFRs = new_dFs2(np.real(np.fft.ifftn(np.fft.ifftshift(e1))), xyz, vol)
    df_val, ddf_val = fcodes_fast.calc_derivatives(
        e0, e1, w_grid, w2_grid, sv_np, dFRs, dRdq, xyz_sum, vol, nx, ny, nz
    )
    ddf_val_inv = np.linalg.pinv(ddf_val)
    step = ddf_val_inv.dot(-df_val)
    # get axis derivatives
    df_ax = np.zeros(3, dtype="float")
    ddf_ax = np.zeros((3, 3), dtype="float")
    dqda = get_dqda(q)
    for i in range(len(dqda)):
        df_ax[i] = df_val[3 + i] * dqda[i]
        for j in range(len(dqda)):
            ddf_ax[i, j] = ddf_val[3 + i, 3 + j] * dqda[i] * dqda[j]
    ddf_ax_inv = np.linalg.pinv(ddf_ax)
    step_ax = ddf_ax_inv.dot(-df_ax)
    return step[:3], step_ax


def derivatives_translation(e0, e1, wgrid, w2grid, sv):
    PI = np.pi
    tp2 = (2.0 * PI) ** 2
    tpi = 2.0 * PI * 1j
    # translation derivatives
    df = np.zeros(3, dtype="float")
    ddf = np.zeros((3, 3), dtype="float")
    for i in range(3):
        df[i] = np.real(
            np.sum(wgrid * e0 * np.conjugate(e1 * tpi * sv[i, :, :, :])))
        for j in range(3):
            if i == 0 or (i > 0 and j >= i):
                ddf[i, j] = -tp2 * \
                    np.sum(w2grid * sv[i, :, :, :] * sv[j, :, :, :])
            else:
                ddf[i, j] = ddf[j, i]
    ddf_inv = np.linalg.pinv(ddf)
    step = ddf_inv.dot(-df)
    return step


def getaxis(axis_ini, step_ax):
    axmax = np.argmax(abs(axis_ini))
    if axmax == 0:
        ay, az = axis_ini[1] + step_ax[1], axis_ini[2] + step_ax[2]
        ax = np.sqrt(1.0 - ay * ay - az * az)
    elif axmax == 1:
        ax, az = axis_ini[0] + step_ax[0], axis_ini[2] + step_ax[2]
        ay = np.sqrt(1.0 - ax * ax - az * az)
    elif axmax == 2:
        ax, ay = axis_ini[0] + step_ax[0], axis_ini[1] + step_ax[1]
        az = np.sqrt(1.0 - ax * ax - ay * ay)
    return np.array([ax, ay, az], dtype="float")


class linefit2:
    def __init__(self):
        self.e0 = None
        self.e1 = None
        self.bin_idx = None
        self.nbin = None
        self.step = None
        self.q_prev = None
        self.t = None
        self.q_init = np.array([1.0, 0.0, 0.0, 0.0], dtype=np.float64)

    def get_linefit_static_data(self, e_list, bin_idx, res_arr, smax):
        if len(e_list) == 2:
            eout, self.bin_idx, self.nbin = cut_resolution_for_linefit(
                e_list, bin_idx, res_arr, smax
            )
        else:
            print("len(e_list: ", len(e_list))
            raise SystemExit()
        self.e0 = eout[0, :, :, :]
        self.e1 = eout[1, :, :, :]

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
        t = self.step * i
        st, _, _, _ = fcodes_fast.get_st(nx, ny, nz, t)
        e1_t = self.e1 * st
        w_grid = self.get_fsc_wght(self.e0, e1_t, self.bin_idx, self.nbin)
        fval = np.sum(w_grid * self.e0 * np.conjugate(e1_t))
        return -fval.real

    def scalar_opt_trans(self):
        f = self.func_t
        res = minimize_scalar(f, method="brent")
        return res.x


class EmFit:
    def __init__(self, mapobj, interp="linear", dfs=None):
        self.mapobj = mapobj
        self.ful_dim = mapobj.map_dim
        self.cell = mapobj.map_unit_cell
        self.pixsize = mapobj.pixsize
        self.origin = mapobj.map_origin
        self.interp = interp
        self.dfs = dfs
        self.cut_dim = mapobj.cdim  # mapobj.map_dim
        self.nbin = mapobj.cbin  # mapobj.nbin
        self.bin_idx = mapobj.cbin_idx  # mapobj.bin_idx
        self.res_arr = mapobj.res_arr[: self.nbin]  # mapobj.res_arr
        self.w_grid = None
        self.fsc = None
        self.sv = None
        self.t = None
        self.st = None
        self.step = None
        self.q = None
        self.axis = None
        self.angle = None
        self.rotmat = None
        self.t_accum = None
        self.ert = None
        self.frt = None
        self.e0 = None
        self.e1 = None
        self.w2_grid = None
        self.fsc_lst = []
        self.avg_fsc = None
        self.q_final_list = []
        self.t_final_list = []
        self.axis_final_list = []

    def calc_fsc(self, e1):
        fsc = core.fsc.anytwomaps_fsc_covariance(self.e0, e1, self.bin_idx, self.nbin)[
            0
        ]
        return fsc

    def calc_fsc_t(self, fstatic, frotated, t):
        cx, cy, cz = fstatic.shape
        st, s1, s2, s3 = fcodes_fast.get_st(cx, cy, cz, t)
        sv = np.array([s1, s2, s3])
        frt = frotated * st
        fsc = core.fsc.anytwomaps_fsc_covariance(fstatic, frt, self.bin_idx, self.nbin)[
            0
        ]
        return frt, fsc, sv

    def rotatemap(self, rotmat, t, smax_lf):
        estatic = self.mapobj.ceo_lst[0]
        fstatic = self.mapobj.cfo_lst[0]
        erotated = get_FRS(rotmat, estatic, interp=self.interp)[:, :, :, 0]
        frotated = get_FRS(rotmat, fstatic, interp=self.interp)[:, :, :, 0]
        # optimize translation
        for i in range(1):
            ert, fsc, sv = self.calc_fsc_t(estatic, erotated, t)
            w_grid, w2_grid = self.get_wght(fsc)
            fval = np.sum(w_grid * estatic * np.conjugate(ert))
            # print(i, fval.real, t)
            step = derivatives_translation(estatic, ert, w_grid, w2_grid, sv)
            lft = linefit2()
            lft.get_linefit_static_data(
                [estatic, ert], self.bin_idx, self.mapobj.res_arr, smax_lf
            )
            lft.step = step
            alpha = lft.scalar_opt_trans()
            t = t + step * alpha
        cx, cy, cz = fstatic.shape
        st, s1, s2, s3 = fcodes_fast.get_st(cx, cy, cz, t)
        ert = erotated * st
        return ert, sv, t

    def get_wght(self, fsc):
        cx, cy, cz = self.e0.shape
        w_grid = fcodes_fast.read_into_grid(
            self.bin_idx,
            fsc,  # fsc / (1 - fsc ** 2),
            self.nbin,
            cx,
            cy,
            cz,
        )
        fsc_sqd = fsc ** 2
        fsc_combi = fsc_sqd,  # fsc_sqd / (1 - fsc_sqd)
        w2_grid = fcodes_fast.read_into_grid(
            self.bin_idx, fsc_combi, self.nbin, cx, cy, cz
        )
        return w_grid, w2_grid

    def functional(self, w_grid, ert):
        fval = np.sum(w_grid * self.e0 * np.conjugate(ert))
        return fval.real

    def get_quaternion(self, angle_i):
        rv = self.axis
        s = math.sin(angle_i / 2.0)
        q1, q2, q3 = rv[0] * s, rv[1] * s, rv[2] * s
        q0 = math.cos(angle_i / 2.0)
        q = np.array([q0, q1, q2, q3], dtype=np.float64)
        return q

    def minimizer(self, ncycles, t_init, q_init_list, smax_lf):
        fsc_lst = []
        fval_list = []
        q_list = []
        t_list = []
        nfit = len(q_init_list)
        print("Refinement started...")
        self.e0 = self.mapobj.ceo_lst[0]
        self.e1 = self.e0
        fstatic = self.mapobj.cfo_lst[0]
        if mapout:
            print("Saving static map")
            data2write = np.real(ifftshift(ifftn(ifftshift(fstatic))))
            core.iotools.write_mrc(data2write, "static_map.mrc", self.cell)
        xyz = create_xyz_grid(self.cell, self.cut_dim)
        vol = self.cut_dim[0] * self.cut_dim[0] * self.cut_dim[0]
        xyz_sum = get_xyz_sum(xyz)
        print("Cycle#   ", "Fval  ", "avg(FSC)")
        self.t = np.asarray(t_init, dtype="float")
        angle_list = []
        nfit = 1
        q = q_init_list[0]
        axis_ang = np.asarray(quaternions.quart2axis(q), dtype="float")
        angle_list.append(axis_ang[-1])
        print("Initial angle: ", axis_ang[-1])
        axis_ini = np.array(axis_ang[:3], dtype="float")
        self.axis = axis_ini
        rotmat = quaternions.get_RM(q)
        for i in range(ncycles):
            e1, sv, self.t = self.rotatemap(rotmat, self.t, smax_lf)
            fsc = self.calc_fsc(e1)
            avg_fsc = np.average(fsc)
            # print(fsc)
            w_grid, w2_grid = self.get_wght(fsc)
            fval = self.functional(w_grid, e1)
            fval_list.append(fval)
            print("{:5d} {:8.4f} {:6.4f}".format(i, fval, avg_fsc))
            print(self.axis)
            if i > 1 and avg_fsc < 0.2:
                print(
                    "Average FSC is lower than 0.4. The axis being refined may not be correct"
                )
                break
            _, step_ax = derivatives(
                self.e0,
                e1,
                w_grid,
                w2_grid,
                q,
                sv,
                xyz,
                xyz_sum,
                vol,
            )
            lft = linefit()
            lft.angle_list = angle_list
            lft.axis_ini = self.axis
            lft.step_axis = step_ax
            lft.get_linefit_static_data(
                [self.e0, self.e0], self.bin_idx, self.res_arr, smax_lf
            )
            alpha = lft.scalar_opt()
            step_ax = step_ax * alpha
            print(step_ax)
            self.axis = getaxis(self.axis, step_ax)
            q = self.get_quaternion(angle_list[0])
            rotmat = quaternions.get_RM(q)
        # Final FSC plot
        print("Refinement finished.")
        q = self.get_quaternion(angle_list[0])
        rotmat = quaternions.get_RM(q)
        cx, cy, cz = fstatic.shape
        st, _, _, _ = fcodes_fast.get_st(cx, cy, cz, self.t)
        frotated = get_FRS(rotmat, fstatic, interp=self.interp)[:, :, :, 0]
        if mapout:
            data2write = np.real(ifftshift(ifftn(ifftshift(frotated * st))))
            outmap = emdcode + "_symrefined_map_fold" + str(fold) + ".mrc"
            core.iotools.write_mrc(data2write, outmap, self.cell)
            print("Written ", outmap)
        print("Now calculating final FSC...")
        fsc = core.fsc.anytwomaps_fsc_covariance(
            fstatic, frotated * st, self.bin_idx, self.nbin
        )[0]
        self.avg_fsc = np.average(fsc)
        fscfile = emdcode + "_fsc" + str(fold) + ".eps"
        core.plotter.plot_nlines(
            self.res_arr,
            [fsc],
            fscfile,
            ["FSC after"],
        )
        print("FSC file was written: ", fscfile)
        # plot fval vs ncycle for clarity
        fig = plt.figure(figsize=(6, 4))
        ax1 = fig.add_subplot(111)
        ax1.plot(fval_list, label="fval", linewidth=2)
        plt.savefig("fval_vs_cycle.eps", format="eps", dpi=300)


def fsc_between_static_and_transfomed_map(
    staticmap, movingmap, bin_idx, rm, t, cell, nbin
):
    nx, ny, nz = staticmap.shape
    st, _, _, _ = fcodes_fast.get_st(nx, ny, nz, t)
    frt_full = get_FRS(rm, movingmap * st, interp="linear")[:, :, :, 0]
    f1f2_fsc = core.fsc.anytwomaps_fsc_covariance(
        staticmap, frt_full, bin_idx, nbin)[0]
    return f1f2_fsc


def get_ibin(bin_fsc, thresh=0.5):
    print('thresh:', thresh)
    # new search from rear end
    for i, ifsc in reversed(list(enumerate(bin_fsc))):
        if ifsc > thresh:
            ibin = i
            if ibin % 2 != 0:
                ibin = ibin - 1
            break
    return ibin


# frequency marching
def run_fit(
    emmap1,
    rotmat,
    t,
    ncycles,
    ifit,
    fitfsc=0.5,
    nmarchingcycles=10,
    fobj=None,
    interp=None,
    fitres=None,
):

    q_init = quaternions.rot2quart(rotmat)
    axis_ang = quaternions.quart2axis(q_init)
    axis_ini = axis_ang[:3]
    angle = axis_ang[-1]
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
    for i in range(nmarchingcycles):
        if i == 0:
            f1f2_fsc = fsc_between_static_and_transfomed_map(
                staticmap=emmap1.fo_lst[0],
                movingmap=emmap1.fo_lst[ifit],
                bin_idx=emmap1.bin_idx,
                rm=rotmat,
                t=t,
                cell=emmap1.map_unit_cell,
                nbin=emmap1.nbin,
            )
            fsc_lst.append(f1f2_fsc)
            # if np.average(f1f2_fsc) > 0.999:
            if fitfsc > 0.999:
                rotmat = rotmat
                final_tran = t
                final_axis = axis_ini
                avg_fsc = fitfsc  # np.average(f1f2_fsc)
                print("\n***FSC between static and moving maps***\n")
                print("bin#     resolution(A)      start-FSC     end-FSC\n")
                for j in range(len(emmap1.res_arr)):
                    print(
                        "{:5d} {:6.2f} {:8.4f} {:8.4f}".format(
                            j, emmap1.res_arr[j], fsc_lst[0][j], fsc_lst[0][j]
                        )
                    )
                break
            ibin = get_ibin(f1f2_fsc, thresh=fitfsc)
            if fitbin < ibin:
                ibin = fitbin
            ibin_old = ibin
            q = q_init
            print("Fitting starts at ", emmap1.res_arr[ibin], " (A)")
        else:
            # Apply initial rotation and translation to calculate fsc
            f1f2_fsc = fsc_between_static_and_transfomed_map(
                emmap1.fo_lst[0],
                emmap1.fo_lst[ifit],
                emmap1.bin_idx,
                rotmat,
                t,
                emmap1.map_unit_cell,
                emmap1.nbin,
            )
            ibin = get_ibin(f1f2_fsc, thresh=fitfsc)
            if fitbin < ibin:
                ibin = fitbin
            if ibin_old == ibin:
                fsc_lst.append(f1f2_fsc)
                q_final = quaternions.rot2quart(rotmat)
                res_arr = emmap1.res_arr[:ibin_old]
                fsc_bef = fsc_lst[0][:ibin_old]
                fsc_aft = fsc_lst[1][:ibin_old]
                print("\n***FSC between static and moving maps***\n")
                print("bin#     resolution(A)      start-FSC     end-FSC\n")
                for j in range(len(res_arr)):
                    print(
                        "{:5d} {:6.2f} {:8.4f} {:8.4f}".format(
                            j, res_arr[j], fsc_bef[j], fsc_aft[j]
                        )
                    )
                break
            else:
                ibin_old = ibin
                print("Fitting starts at ", emmap1.res_arr[ibin], " (A)")
        if ibin == 0:
            print("ibin = 0, Cannot proceed! Stopping current axis refinement.")
            fobj.write(
                "ibin = 0, Cannot proceed! Stopping current axis refinement.\n")
            break
        e_list = [emmap1.eo_lst[0], emmap1.fo_lst[0]]
        eout, cBIdx, cbin = cut_resolution_for_linefit(
            e_list, emmap1.bin_idx, emmap1.res_arr, ibin
        )
        #static_cutmap = eout[0, :, :, :]
        static_cutmap = eout[1, :, :, :]  # use Fo instead of Eo for fitting.
        fstatic_cutmap = eout[1, :, :, :]
        moving_cutmap = static_cutmap
        emmap1.ceo_lst = [static_cutmap, moving_cutmap]
        emmap1.cfo_lst = [fstatic_cutmap]
        emmap1.cbin_idx = cBIdx
        emmap1.cdim = moving_cutmap.shape
        emmap1.cbin = cbin
        fit = EmFit(emmap1)
        slf = ibin
        fit.minimizer(ncycles=ncycles, t_init=t, smax_lf=slf, q_init_list=[q])
        ncycles = ncycles
        final_axis = fit.axis
        final_tran = fit.t
        avg_fsc = fit.avg_fsc
        q = quaternions.get_quaternion(
            [[final_axis[0], final_axis[1], final_axis[2]],
                [np.rad2deg(angle)]]
        )
        rotmat = quaternions.get_RM(q)
    return final_axis, final_tran, avg_fsc


def overlay(
    emdcode,
    imap,
    rotaxis,
    symorder,
    fitfsc=0.5,
    ncycles=10,
    t_init=[0.0, 0.0, 0.0],
    theta_init=0.0,
    smax=7,
    interp="linear",
    imask=None,
    fobj=None,
    halfmaps=False,
    dfs_interp=False,
    usemodel=False,
    fitres=5,
):
    axis = np.asarray(rotaxis)
    axis = axis / math.sqrt(np.dot(axis, axis))
    initial_axis = axis
    print("Initial axis and fold: ", axis, symorder)
    fobj.write("Initial axis and fold: " + str(axis) + str(symorder) + "\n")
    print(" ")
    print("Preparing data for axis refinement...")
    try:
        emmap1 = EmmapOverlay(imap=imap, imask=imask)
    except:
        emmap1 = EmmapOverlay(imap=imap)
    emmap1.get_maps()
    emmap1.calc_fsc_from_maps()
    print(" Number of refinement cycles:", ncycles)
    print("Data resolution for refinement: ", fitres)
    fobj.write(" Number of refinement cycles: " + str(ncycles) + "\n")
    fobj.write("Data resolution for refinement: " + str(fitres) + "\n")
    q_init_list = []
    q_final_list = []
    avg_fsc_list = []
    angle_list = []
    rotmat_list = []
    axis_list = []
    trans_list = []
    print("Initial axis and angles:")
    fobj.write("Initial axis and angles: \n")
    for n in range(1, symorder):
        theta = float(n * 360.0 / symorder)
        print("   ", axis, theta)
        fobj.write("   " + str(axis) + str(theta) + "\n")
        angle_list.append(theta)
        q = quaternions.get_quaternion([[axis[2], axis[1], axis[0]], [theta]])
        # print(q)
        q_init_list.append(q)

    for ifit in range(len(emmap1.eo_lst)):
        rotmat_init = quaternions.get_RM(q_init_list[ifit])
        final_axis, final_tran, avg_fsc = run_fit(
            emmap1=emmap1,
            rotmat=rotmat_init,
            t=np.asarray(t_init, dtype="float"),
            ncycles=ncycles,
            ifit=ifit,
            interp=interp,
            fitres=fitres,
            fobj=fobj,
            fitfsc=fitfsc,
        )
        axis_list.append(final_axis)
        trans_list.append(final_tran)
        avg_fsc_list.append(avg_fsc)
    for ax, tr, avgfsc in zip(axis_list, trans_list, avg_fsc_list):
        final_axis = ax
        final_tran = tr
        # print("Final axis: ", final_axis[::-1])
        # print("Final translation (A): ", final_tran * emmap1.pixsize)
        # fobj.write("   Final axis: " + str(final_axis[::-1]) + "\n")
    return initial_axis, final_axis[::-1], avgfsc


def get_intial_axis(imap):
    #fold, x, y, z, peakh = proshade.get_symmops_from_proshade(imap)
    results = proshade.get_symmops_from_proshade(imap)
    return results


def prime_factors(n):
    # https://stackoverflow.com/questions/15347174/python-finding-prime-factors
    i = 2
    factors = []
    while i * i <= n:
        if n % i:
            i += 1
        else:
            n //= i
            factors.append(i)
    if n > 1:
        factors.append(n)
    return factors


def filter_axes(imap, resol, use_proshade_peakheight=True, use_fsc=False,
                peak_cutoff=0.8, fsc_cutoff=0.7, ang_tol=5.0, fobj=None):

    uc, arr, orig = em.get_data(imap)
    #arr, uc = reboxmap(arr=arr, uc=uc)
    arr = set_dim_even(arr)
    f1 = fftshift(fftn(fftshift(arr)))
    nbin, res_arr, bin_idx = core.restools.get_resolution_array(uc, f1)
    resol = float(resol) * 1.1
    dist = np.sqrt((res_arr - resol) ** 2)
    ibin = np.argmin(dist)
    # get initial axes list from ProSHADE
    proshade_results = get_intial_axis(imap)
    symorder = proshade_results[0]
    x = proshade_results[1]
    y = proshade_results[2]
    z = proshade_results[3]
    peakh = proshade_results[4]
    proshade_pg = proshade_results[5]
    #symorder, x, y, z, peakh = get_intial_axis(imap)
    fobj.write("ProSHADE peak table \n")
    if len(symorder) < 1:
        print("proshade peak table is empty.")
        return []
        #SystemExit()
    for i, odr in enumerate(symorder):
        fobj.write(str(symorder[i]) +" ["+ str(x[i]) +" "+ str(y[i]) 
            +" "+ str(z[i]) +"] "+ str(peakh[i]) + "\n")
    fobj.write("Proshade Point group: "+ proshade_pg + "\n")
    #
    # find the peak_cutoff if igen cutoff is too high
    peakcutoff_np = np.asarray(peakh, 'float')
    peakcutoff_best = np.max(peakcutoff_np) - 0.1
    if peak_cutoff > np.max(peakcutoff_np):
        peak_cutoff = peakcutoff_best
        print("new peak cutoff: ", peak_cutoff)
        fobj.write("new peak cutoff: "+ str(peak_cutoff) + "\n")
    axes_list = []
    order_list = []
    avgfsc_list = []
    if use_proshade_peakheight and use_fsc:
        use_proshade_peakheight = False
    if use_proshade_peakheight:
        use_fsc = False
        for ix, iy, iz, order, score in zip(x, y, z, symorder, peakh):
            # using proshade peak height to decide point group
            if score > peak_cutoff:
                axis = np.asarray([ix, iy, iz], 'float')
                axis = axis / math.sqrt(np.dot(axis, axis))
                axes_list.append(axis)
                order_list.append(order)
                avgfsc_list.append(score)
    if use_fsc:
        print("ORDER     FSC(avg)       P.H      AXIS")
        fobj.write("ORDER     FSC(avg)       P.H      AXIS \n")
        for ix, iy, iz, order, score in zip(x, y, z, symorder, peakh):
            # using proshade peak height to decide point group
            if score > 0.1:
                axis = np.asarray([ix, iy, iz], 'float')
                axis = axis / math.sqrt(np.dot(axis, axis))
                # calculate FSC, because score is not good enough to decide the pg
                theta = float(360.0 / order)
                q = quaternions.get_quaternion(
                    [[axis[2], axis[1], axis[0]], [theta]])
                rotmat = quaternions.get_RM(q)
                frt = get_FRS(rotmat, f1, interp="linear")[:, :, :, 0]
                fsc = core.fsc.anytwomaps_fsc_covariance(
                    f1, frt, bin_idx, nbin)[0]
                avg_fsc = np.average(fsc[:ibin])
                print(order, avg_fsc, score, axis)
                fobj.write(str(order) +" "+ str(avg_fsc) +" "+ str(score) +" "+ str(axis) + "\n")
                if avg_fsc > fsc_cutoff:
                    axes_list.append(axis)
                    order_list.append(order)
                    avgfsc_list.append(avg_fsc)
    # sort all lists by order
    sorted_list = sort_together([order_list, axes_list, avgfsc_list])
    sorted_odrlist, sorted_axlist, sorted_fsclist = sorted_list
    # remove duplicate axes
    duplicate_list = []
    for i, odr1 in enumerate(sorted_odrlist):
        for j, odr2 in enumerate(sorted_odrlist):
            if i < j:
                angle = cosine_angle(sorted_axlist[i], sorted_axlist[j])
                if angle < ang_tol or (180. - angle) < ang_tol:
                    if odr1 == odr2:
                        duplicate_list.append(j)
                    elif odr1 < odr2:
                        duplicate_list.append(i)
                    elif odr1 > odr2:
                        duplicate_list.append(j)
    print(duplicate_list)
    cleaned_axlist = []
    cleaned_odrlist = []
    cleaned_fsclist = []
    for i, ax in enumerate(sorted_axlist):
        if not np.any(i == np.asarray(duplicate_list)):
            print(i, ax, sorted_odrlist[i], sorted_fsclist[i])
            cleaned_axlist.append(ax)
            cleaned_odrlist.append(sorted_odrlist[i])
            cleaned_fsclist.append(sorted_fsclist[i])
    return [cleaned_axlist, cleaned_odrlist, cleaned_fsclist, f1, ibin, nbin, bin_idx, proshade_pg]


def prefilter_order(f1, axis, order, ibin, nbin, bin_idx):
    # prime factorization
    factors = prime_factors(order)
    filtered_order_list = []
    if len(factors) == 1:
        filtered_order_list.append(order)
    else:
        axis = np.asarray(axis)
        axis = axis / math.sqrt(np.dot(axis, axis))
        rotmat_list = []
        for fac in factors:
            theta = float(360.0 / fac)
            q = quaternions.get_quaternion(
                [[axis[2], axis[1], axis[0]], [theta]])
            rotmat = quaternions.get_RM(q)
            frt = get_FRS(rotmat, f1, interp="linear")[:, :, :, 0]
            fsc = core.fsc.anytwomaps_fsc_covariance(f1, frt, bin_idx, nbin)[0]
            avg_fsc = np.average(fsc[:ibin])
            print('fac, Avg FSC: ', fac, avg_fsc)
            if avg_fsc > 0.8:
                filtered_order_list.append(fac)
    true_order = np.prod(np.asarray(filtered_order_list), dtype='int')
    return true_order


""" def prefilter_order(imap, axis, order, resol=None):
    # prime factorization
    factors = prime_factors(order)
    filtered_order_list = []
    if len(factors) == 1:
        filtered_order_list.append(order)
    else:
        axis = np.asarray(axis)
        axis = axis / math.sqrt(np.dot(axis, axis))
        uc, arr, orig = em.get_data(imap)
        arr = set_dim_even(arr)
        f1 = fftshift(fftn(fftshift(arr)))
        nbin, res_arr, bin_idx = core.restools.get_resolution_array(uc, f1)
        rotmat_list = []
        if resol is not None:
            dist = np.sqrt((res_arr - resol) ** 2)
            ibin = np.argmin(dist)
            # cut map
            #f1, bin_idx, nbin = cut_resolution_for_linefit(
            #[f1], bin_idx, res_arr, ibin)
        else:
            ibin = nbin

        for fac in factors:
            theta = float(360.0 / fac)
            q = quaternions.get_quaternion([[axis[2], axis[1], axis[0]], [theta]])
            rotmat = quaternions.get_RM(q)
            frt = get_FRS(rotmat, f1, interp="linear")[:, :, :, 0]
            fsc = core.fsc.anytwomaps_fsc_covariance(f1, frt, bin_idx, nbin)[0]
            avg_fsc = np.average(fsc[:ibin])
            print('fac, Avg FSC: ', fac, avg_fsc)
            if avg_fsc > 0.8:
                filtered_order_list.append(fac)
    true_order = np.prod(np.asarray(filtered_order_list), dtype='int')
    return true_order """


def cosine_angle(ax1, ax2):
    vec_a = np.asarray(ax1, 'float')
    vec_b = np.asarray(ax2, 'float')
    dotp = np.dot(vec_a, vec_b)
    #assert -1.0 <= dotp <= 1.0
    if -1.0 <= dotp <= 1.0:
        angle = math.acos(dotp)
    else:
        print("Problem, dotp: ", dotp)
        print("axes:", vec_a, vec_b)
        angle = 0.0
    return np.rad2deg(angle)


def decide_pointgroup(axeslist, orderlist):
    # TODO- Current group generator axes are not correct.
    # need to identify correct axes for refinement.

    # check for cyclic sym of n-order
    order_arr = np.asarray(orderlist)
    dic = {i: (order_arr == i).nonzero()[0] for i in np.unique(order_arr)}
    uniqorder = np.fromiter(dic.keys(), dtype='int')
    #uniqorder = np.unique(np.asarray(orderlist), 'int')
    anglist = []
    point_group = None
    gp_generator_ax1 = None
    gp_generator_ax2 = None
    order1 = order2 = 0
    ang_tol = 5.0  # Degrees
    print("uniqorders: ", uniqorder)
    if len(uniqorder) == 1:
        print("Len of uniqorder: ", 1)
        if uniqorder[0] == 1:
            print("Unsymmetrized map")
            point_group = 'C1'
            order1 = 1
        elif uniqorder[0] != 1:
            odrn = dic[uniqorder[0]]
            if len(odrn) == 1:
                point_group = 'C'+str(uniqorder[0])
                gp_generator_ax1 = axeslist[0]
                order1 = uniqorder[0]
            elif len(odrn) > 1:
                if uniqorder[0] == 2:
                    print("Checking for D symmetry...")
                    for i in odrn:
                        for j in odrn:
                            if i != j:
                                angle = cosine_angle(axeslist[i], axeslist[j])
                                anglist.append(angle)
                    angarr = np.asarray(anglist, 'float')
                    if np.all(abs(angarr - 90.0) <= ang_tol):
                        point_group = 'D2'
                        gp_generator_ax1 = axeslist[odrn[0]]
                        gp_generator_ax2 = axeslist[odrn[1]]
                        order1 = order2 = 2
                    else:
                        print("Unknown symmetry")
                        point_group = 'Unkown'
                elif uniqorder[0] == 3:
                    print("Checking for T symmetry...")
                    for i in odrn:
                        for j in odrn:
                            if i != j:
                                angle = cosine_angle(axeslist[i], axeslist[j])
                                anglist.append(angle)
                    angarr = np.asarray(anglist, 'float')
                    condition2 = angarr[np.logical_not(
                        abs(angarr - 109.47) <= ang_tol)]
                    if np.all(abs(condition2 - 70.53) <= ang_tol):
                        point_group = 'T'
                        gp_generator_ax1 = axeslist[odrn[0]]
                        gp_generator_ax2 = axeslist[odrn[1]]
                        order1 = order2 = 3
                    else:
                        print("test1")
                        print("Unknown point group.")
                        point_group = 'Unkown'
                elif uniqorder[0] == 4:
                    print("Checking for O symmetry...")
                    for i in odrn:
                        for j in odrn:
                            if i != j:
                                angle = cosine_angle(axeslist[i], axeslist[j])
                                anglist.append(angle)
                    angarr = np.asarray(anglist, 'float')
                    if np.all(abs(angarr - 90.0) <= ang_tol):
                        point_group = 'O'
                        gp_generator_ax1 = axeslist[odrn[0]]
                        gp_generator_ax2 = axeslist[odrn[1]]
                        order1 = order2 = 4
                    else:
                        print("Unknown symmetry")
                        point_group = 'Unkown'
                elif uniqorder[0] == 5:
                    print("Checking for I symmetry...")
                    for i in odrn:
                        for j in odrn:
                            if i != j:
                                angle = cosine_angle(axeslist[i], axeslist[j])
                                anglist.append(angle)
                    angarr = np.asarray(anglist, 'float')
                    condition2 = angarr[np.logical_not(
                        abs(angarr - 63.47) <= ang_tol)]
                    if np.all(abs(condition2 - 116.57) <= ang_tol):
                        point_group = 'I'
                        gp_generator_ax1 = axeslist[odrn[0]]
                        gp_generator_ax2 = axeslist[odrn[1]]
                        order1 = order2 = 5
                    else:
                        print("Unknown symmetry")
                        point_group = 'Unkown'
                else:
                    print("Unknown symmetry")
                    point_group = 'Unkown'
    elif len(uniqorder) == 2:
        if np.any(uniqorder == 2):
            odr2 = dic[2]
            if np.any(uniqorder == 3):
                odr3 = dic[3]  # get all 3-fold axes locations
                if len(odr3) == 1:
                    print("Ckecking for D symmetry...")
                    for i in odr2:
                        for j in odr3:
                            angle = cosine_angle(axeslist[i], axeslist[j])
                            anglist.append(angle)
                    angarr = np.asarray(anglist, 'float')
                    if np.all(abs(angarr - 90.0) <= ang_tol):
                        point_group = 'D3'
                        gp_generator_ax1 = axeslist[odr2[0]]
                        gp_generator_ax2 = axeslist[odr3[0]]
                        order1 = 2
                        order2 = 3
                    else:
                        print("Unknown symmetry")
                        point_group = 'Unkown'
                if len(odr3) > 1:
                    print("Checking for T symmetry...")
                    for i in odr3:
                        for j in odr3:
                            if i != j:
                                angle = cosine_angle(axeslist[i], axeslist[j])
                                anglist.append(angle)
                    angarr = np.asarray(anglist, 'float')
                    #condition1 = angarr[abs(angarr - 109.47) <= ang_tol]
                    condition2 = angarr[np.logical_not(
                        abs(angarr - 109.47) <= ang_tol)]
                    if np.all(abs(condition2 - 70.53) <= ang_tol):
                        point_group = 'T'
                        gp_generator_ax1 = axeslist[odr3[0]]
                        gp_generator_ax2 = axeslist[odr3[1]]
                        order1 = order2 = 3
                    else:
                        print("test2")
                        print("Unknown point group.")
                        point_group = 'Unkown'
            elif np.any(uniqorder == 4):
                odr4 = dic[4]
                if len(odr4) == 1:
                    print("Checking for D symmetry...")
                    for i in odr2:
                        for j in odr4:
                            angle = cosine_angle(axeslist[i], axeslist[j])
                            anglist.append(angle)
                    angarr = np.asarray(anglist, 'float')
                    if np.all(abs(angarr - 90.0) <= ang_tol):
                        point_group = 'D4'
                        gp_generator_ax1 = axeslist[odr2[0]]
                        gp_generator_ax2 = axeslist[odr4[0]]
                        order1 = 2
                        order2 = 4
                    else:
                        print("Unknown symmetry")
                        point_group = 'Unkown'
                elif len(odr4) > 1:
                    print("Checking for O symmetry...")
                    for i in odr4:
                        for j in odr4:
                            if i != j:
                                angle = cosine_angle(axeslist[i], axeslist[j])
                                anglist.append(angle)
                    angarr = np.asarray(anglist, 'float')
                    if np.all(abs(angarr - 90.0) <= ang_tol):
                        point_group = 'O'
                        gp_generator_ax1 = axeslist[odr4[0]]
                        gp_generator_ax2 = axeslist[odr4[1]]
                        order1 = order2 = 4
                    else:
                        print("Unknown symmetry")
                        point_group = 'Unkown'
            elif np.any(uniqorder == 5):
                odr5 = dic[5]
                if len(odr5) == 1:
                    print("Checking for D symmetry...")
                    for i in odr2:
                        for j in odr5:
                            angle = cosine_angle(axeslist[i], axeslist[j])
                            anglist.append(angle)
                    angarr = np.asarray(anglist, 'float')
                    if np.all(abs(angarr - 90.0) <= ang_tol):
                        point_group = 'D5'
                        gp_generator_ax1 = axeslist[odr2[0]]
                        gp_generator_ax2 = axeslist[odr5[0]]
                        order1 = 2
                        order2 = 5
                    else:
                        print("Unknown symmetry")
                        point_group = 'Unkown'
                if len(odr5) > 1:
                    print("Checking for I symmetry...")
                    for i in odr5:
                        for j in odr5:
                            if i != j:
                                angle = cosine_angle(axeslist[i], axeslist[j])
                                anglist.append(angle)
                    angarr = np.asarray(anglist, 'float')
                    condition2 = angarr[np.logical_not(
                        abs(angarr - 63.47) <= ang_tol)]
                    if np.all(abs(condition2 - 116.57) <= ang_tol):
                        point_group = 'I'
                        gp_generator_ax1 = axeslist[odr5[0]]
                        gp_generator_ax2 = axeslist[odr5[1]]
                        order1 = order2 = 5
                    else:
                        print("Unknown symmetry")
                        point_group = 'Unkown'
            else:
                n = uniqorder[uniqorder != 2][0]
                odrn = dic[n]
                print("Ckecking for D symmetry...")
                for i in odr2:
                    for j in odrn:
                        angle = cosine_angle(axeslist[i], axeslist[j])
                        anglist.append(angle)
                angarr = np.asarray(anglist, 'float')
                if np.all(abs(angarr - 90.0) <= ang_tol):
                    point_group = 'D' + str(n)
                    gp_generator_ax1 = axeslist[odr2[0]]
                    gp_generator_ax2 = axeslist[odrn[0]]
                    order1 = 2
                    order2 = n
                else:
                    print("Unknown symmetry")
                    point_group = 'Unkown'
    elif len(uniqorder) > 2:
        # groups must belong to O or I
        if np.any(uniqorder == 2):
            odr2 = dic[2]
            if np.any(uniqorder == 3):
                odr3 = dic[3]
                if np.any(uniqorder == 4):
                    odr4 = dic[4]
                    print("Ckecking for O symmetry...")
                    for i in odr4:
                        for j in odr4:
                            if i != j:
                                angle = cosine_angle(axeslist[i], axeslist[j])
                                anglist.append(angle)
                    angarr = np.asarray(anglist, 'float')
                    if np.all(abs(angarr - 90.0) <= ang_tol):
                        point_group = 'O'
                        gp_generator_ax1 = axeslist[odr4[0]]
                        gp_generator_ax2 = axeslist[odr4[1]]
                        order1 = order2 = 4
                    else:
                        print("Unknown point group.")
                        point_group = 'Unkown'
                elif np.any(uniqorder == 5):
                    # I symmetry
                    odr5 = dic[5]
                    print("Ckecking for I symmetry...")
                    for i in odr5:
                        for j in odr5:
                            if i != j:
                                angle = cosine_angle(axeslist[i], axeslist[j])
                                anglist.append(angle)
                    angarr = np.asarray(anglist, 'float')
                    #condition1 = angarr[abs(angarr - 63.47) <= ang_tol]
                    condition2 = angarr[np.logical_not(
                        abs(angarr - 63.47) <= ang_tol)]
                    if np.all(abs(condition2 - 116.57) <= ang_tol):
                        point_group = 'I'
                        gp_generator_ax1 = axeslist[odr5[0]]
                        gp_generator_ax2 = axeslist[odr5[1]]
                        order1 = order2 = 5
                    else:
                        print("Unknown point group.")
                        point_group = 'Unkown'
        else:
            print("Unknown symmetry")
            point_group = 'Unkown'
    return point_group, [order1, order2, gp_generator_ax1, gp_generator_ax2]


def get_pg(imap, resol, use_peakheight, peak_cutoff, use_fsc, fsc_cutoff, ang_tol, fobj=None):
    results = filter_axes(imap=imap,
                          resol=resol,
                          use_proshade_peakheight=use_peakheight,
                          use_fsc=use_fsc,
                          peak_cutoff=peak_cutoff,
                          fsc_cutoff=fsc_cutoff,
                          ang_tol=ang_tol,
                          fobj=fobj)
    if len(results) < 1:
        return []
    cleaned_axlist, cleaned_odrlist = results[0], results[1]
    proshade_pg = results[7]
    pg, gp_generators = decide_pointgroup(
        axeslist=cleaned_axlist, orderlist=cleaned_odrlist)
    print("Point group detected from the map: ", pg)
    fobj.write("Point group detected from the map: " + str(pg) + "\n")
    gen_axlist = [gp_generators[2], gp_generators[3]]
    gen_odrlist = [gp_generators[0], gp_generators[1]]
    return [proshade_pg, pg, gen_odrlist, gen_axlist]


def refine_pg_generator_axes(imap, axlist, odrlist, fobj, fitres=5.0, fitfsc=0.7, emdid=None):
    global fold, emdcode
    if emdid is None:
        emdcode = "EMD-0000"
    else:
        emdcode = "EMD-"+str(emdid)
    iniaxlst = []
    fnlaxlst = []
    foldlst = []
    fsclst = []
    for ax, fold in zip(axlist, odrlist):
        if fold != 0:
            foldlst.append(fold)
            initial_axis, final_axis, avg_fsc = overlay(
                emdcode=emdcode, imap=imap, rotaxis=ax, symorder=fold, fobj=fobj, fitfsc=fitfsc, fitres=fitres)
            iniaxlst.append(ax)
            fnlaxlst.append(final_axis)
            fsclst.append(avg_fsc)
    print("**Refined results**")
    print("Initial ax, Refined ax, FSC_avg")
    for i, _ in enumerate(foldlst):
        print(foldlst[i], iniaxlst[i], fnlaxlst[i], fsclst[i])
    return [foldlst, fnlaxlst]


def analyse_nmaps(maplist, reslist, use_peakheight=True, peak_cutoff=0.8,
                  use_fsc=False, fsc_cutoff=0.7, ang_tol=5.0, refine_axes=False,emdlist=None):
    fobj = open("EMDA_symref.txt", "w")
    pglist = []
    foldlist = []
    initaxlist = []
    fnlaxlist = []

    assert len(maplist) == len(reslist)
    if emdbidlist is not None:
        assert len(emdbidlist) == len(maplist)
    else:
        ids = np.arange(len(maplist))
        emdbidlist = [",".join(item) for item in ids.astype(str)]

    for i, imap in enumerate(maplist):
        emdcode = "EMD-" + emdlist[i]
        resol = float(reslist[i]) * 1.1
        print("Map: ", imap)
        print("emdcode: ", emdcode)
        fobj.write("\n")
        fobj.write("Map: " + imap + "\n")
        fobj.write("emdcode: " + emdcode + "\n")
        print(" ")
        # get the point group in map
        pg, gen_odrs, gen_axes = get_pg(imap=imap,
                                        resol=resol,
                                        use_peakheight=use_peakheight,
                                        peak_cutoff=peak_cutoff,
                                        use_fsc=use_fsc,
                                        fsc_cutoff=fsc_cutoff,
                                        ang_tol=ang_tol)
        pglist.append(pg)
        initaxlist.append(gen_axes)
        # refine gp generator axes
        if refine_axes:
            folds, fnlaxes = refine_pg_generator_axes(
                imap, gen_axlist, gen_odrlist, fobj, fitres=resol, fitfsc=0.7)
            foldlist.append(folds)
            fnlaxlist.append(fnlaxes)
    if refine_axes:
        return pglist, foldlist, initaxlist, fnlaxlist
    else:
        return pglist


if __name__ == "__main__":
    """ maplist = [
        "/Users/ranganaw/MRC/REFMAC/EMD-0882/map/emd_0882.map",
        "/Users/ranganaw/MRC/REFMAC/EMD-0882/emda_test/group_rotation/test_pipeline/emd_0959.map",
        "/Users/ranganaw/MRC/REFMAC/EMD-6952//map/emd_6952.map",
    ] """
    map1, resol = sys.argv[1:]

    maplist = []
    maplist.append(map1)
    reslist = []
    reslist.append(float(resol))
    print(maplist, reslist)
    main(maplist, reslist)
    #filter_axes(map1, resol)
