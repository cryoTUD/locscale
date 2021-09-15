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
from emda.core import fsc as fsctools
from emda.core import quaternions
from emda.ext.mapfit.utils import get_FRS, create_xyz_grid, get_xyz_sum
import emda.emda_methods as em

np.set_printoptions(suppress=True)  # Suppress insignificant values for clarity

timeit = False


class EmFit:
    def __init__(self, mapobj, interp="linear", dfs=None):
        self.mapobj = mapobj
        self.cut_dim = mapobj.cdim
        self.ful_dim = mapobj.map_dim
        self.cell = mapobj.map_unit_cell
        self.pixsize = mapobj.pixsize
        self.origin = mapobj.map_origin
        self.interp = interp
        self.dfs = dfs
        self.w_grid = None
        self.fsc = None
        self.sv = None
        self.t = None
        self.st = None
        self.step = None
        self.q = None
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
        #self.ert = get_FRS(self.rotmat, self.e1 * self.st, interp=self.interp)[:, :, :, 0]
        self.ert = get_FRS(self.rotmat, self.e1, interp=self.interp)[:, :, :, 0]
        self.ert = self.ert  * self.st
        fsc = fsctools.anytwomaps_fsc_covariance(
            self.e0, self.ert, self.mapobj.cbin_idx, self.mapobj.cbin
        )[0]
        return fsc

    def calc_fsc3d(self):
        from emda.ext.fouriersp_local import get_3d_fouriercorrelation
        from emda.core.restools import get_resArr, remove_edge, create_soft_edged_kernel_pxl
        cx, cy, cz = self.e0.shape
        self.st, s1, s2, s3 = fcodes_fast.get_st(cx, cy, cz, self.t)
        self.sv = np.array([s1, s2, s3])
        self.ert = get_FRS(self.rotmat, self.e1 * self.st, interp=self.interp)[:, :, :, 0]
        kern = create_soft_edged_kernel_pxl(r1=5)
        fsc3d = get_3d_fouriercorrelation(hf1=self.e0, hf2=self.ert, kern_sphere=kern)[0]
        nx, ny, nz = self.e0.shape
        fResArr = get_resArr(self.mapobj.map_unit_cell, nx)
        cut_mask = remove_edge(fResArr, fResArr[-1])
        cc_mask = np.zeros(shape=(nx, ny, nz), dtype="int")
        cx, cy, cz = cut_mask.shape
        dx = (nx - cx) // 2
        dy = (ny - cy) // 2
        dz = (nz - cz) // 2
        cc_mask[dx : dx + cx, dy : dy + cy, dz : dz + cz] = cut_mask
        return fsc3d.real * cc_mask

    def get_wght3d(self, fsc3d):
        wgrid3d = fsc3d / (1 - fsc3d**2)
        w2grid3d = fsc3d**2 / (1 - fsc3d**2)
        return wgrid3d, w2grid3d

    def get_wght(self):
        cx, cy, cz = self.e0.shape
        w_grid = fcodes_fast.read_into_grid(
            self.mapobj.cbin_idx, self.fsc / (1 - self.fsc ** 2), self.mapobj.cbin, cx, cy, cz
        )
        fsc_sqd = self.fsc ** 2
        fsc_combi = fsc_sqd / (1 - fsc_sqd)
        w2_grid = fcodes_fast.read_into_grid(
            self.mapobj.cbin_idx, fsc_combi, self.mapobj.cbin, cx, cy, cz
        )
        return w_grid, w2_grid

    def functional(self):
        fval = np.sum(self.w_grid * self.e0 * np.conjugate(self.ert))
        return fval.real

    def minimizer(self, ncycles, t_init, rotmat, smax_lf, fobj=None):
        import math
        from emda.ext.mapfit import (
            rotmat2quart,
            linefit_class,
            interp_derivatives,
            derivatives,
        )
        tol = 1e-2
        fsc_lst = []
        fval_list = []
        q_list = []
        t_list = []
        nfit = len(self.mapobj.ceo_lst) - 1
        self.e0 = self.mapobj.ceo_lst[0]  # Static map e-data for fit
        fobj.write("\n")
        fobj.write("Normalized Structure Factors are used for fitting! \n")
        xyz = create_xyz_grid(self.cell, self.cut_dim)
        #vol = self.cell[0] * self.cell[1] * self.cell[2]
        vol = self.cut_dim[0] * self.cut_dim[0] * self.cut_dim[0]
        xyz_sum = get_xyz_sum(xyz)
        q_init = np.array([1.0, 0.0, 0.0, 0.0], dtype=np.float64)
        fobj.write("\n")
        fobj.write("Cycle#\n")
        fobj.write("Function value\n")
        fobj.write("Rotation(degrees)\n")
        fobj.write("Translation(A)\n")
        fobj.write("\n")
        print("Cycle#   ", "Fval  ", "Rot(deg)  ", "Trans(A)  ", "avg(FSC)")
        for ifit in range(nfit):
            self.e1 = self.mapobj.ceo_lst[ifit + 1]
            for i in range(ncycles):
                start = timer()
                if i == 0:
                    self.t = np.asarray(t_init, dtype='float')
                    t_accum = self.t
                    t_accum_angstrom = trans_in_angstrom(t_accum, self.pixsize, self.ful_dim)
                    translation_vec = np.sqrt(
                        np.sum(t_accum_angstrom * t_accum_angstrom)
                    )
                    self.rotmat = rotmat
                    q = rotmat2quart.rot2quart(self.rotmat)
                    self.q = q
                    q_accum = self.q
                    theta2 = np.arccos((np.trace(self.rotmat) - 1) / 2) * 180.0 / np.pi
                    t_accum_previous = t_accum
                else:
                    rm_accum = quaternions.get_RM(q_accum)
                    theta2 = np.arccos((np.trace(rm_accum) - 1) / 2) * 180.0 / np.pi
                    if theta2 < 0.01:
                        self.rotmat = np.identity(3)
                    else:
                        self.rotmat = quaternions.get_RM(self.q)
                self.fsc = self.calc_fsc()
                #fsc3d = self.calc_fsc3d() 
                #print(self.fsc)
                #em.write_mrc(fsc3d, "fsc3d.mrc", self.cell)
                if np.average(self.fsc) > 0.999:
                    fval = np.sum(self.e0 * np.conjugate(self.ert))
                    print("fval, FSC_avg ", fval.real, np.average(self.fsc))
                    self.rotmat = quaternions.get_RM(q_accum)
                    self.t_accum = t_accum_previous  # final translation
                    self.fsc_lst = fsc_lst
                    translation_vec = trans_in_angstrom(self.t_accum, self.pixsize, self.ful_dim)
                    theta2 = np.arccos((np.trace(self.rotmat) - 1) / 2) * 180.0 / np.pi
                    print("thets2, trans: ", theta2, translation_vec)
                    break
                self.w_grid, self.w2_grid = self.get_wght()
                #self.w_grid, self.w2_grid = self.get_wght3d(fsc3d)
                fval = self.functional()
                fval_list.append(fval)
                q_list.append(q_accum)
                t_list.append(t_accum)
                if math.isnan(theta2):
                    raise SystemExit("Cannot find a solution! Stopping now...")
                if i == 0:
                    fval_previous = fval
                    fsc = self.fsc
                    fsc_lst.append(fsc)
                if i > 0 and fval_previous < fval or i == ncycles - 1:
                    fsc = self.fsc
                # stop fitting if converged
                if i > 0 and abs(fval - fval_previous) <= tol or i == ncycles - 1:
                    self.rotmat = quaternions.get_RM(q_list[-1])
                    self.q = q_list[-1]
                    self.t_accum = t_list[-1]
                    self.rotmat = quaternions.get_RM(q_list[-1])
                    break
                """ if i > 0 and i == ncycles - 1:
                    # search for max fval in the fval_list
                    self.rotmat = quaternions.get_RM(q_list[-1])
                    self.q = q_list[-1]
                    self.t_accum = t_list[-1]
                    self.fsc_lst = fsc_lst
                    translation_vec = trans_in_angstrom(t_accum, self.pixsize, self.ful_dim)
                    theta2 = np.arccos((np.trace(self.rotmat) - 1) / 2) * 180.0 / np.pi
                    print("thets2, trans: ", theta2, translation_vec)
                    break """
                print(
                    "{:5d} {:8.4f} {:6.2f} {:6.2f} {:6.2f}".format(
                        i, fval, theta2, translation_vec, np.average(self.fsc)
                    )
                )
                fobj.write(
                    "{:5d} {:8.4f} {:6.2f} {:6.2f}\n".format(
                        i, fval, theta2, translation_vec
                    )
                )
                if self.dfs is not None:
                    start1 = timer()
                    dFRs = interp_derivatives.interp_derivatives(self.dfs, self.rotmat)
                    end1 = timer()
                    if timeit:
                        print(" time for dFRs calculation: ", end1 - start1)
                elif self.dfs is None:
                    dFRs = None
                t_accum_previous = t_accum
                self.step, _ = derivatives.new_derivatives(
                    self.e0,
                    self.ert,
                    self.w_grid,
                    self.w2_grid,
                    q,
                    self.sv,
                    xyz,
                    xyz_sum,
                    vol,
                    dFRs,
                )
                self.e1 = self.ert
                start_lf = timer()
                if i == 0:
                    linefit = linefit_class.LineFit()
                    linefit.get_linefit_static_data(
                        self.e0, self.mapobj.cbin_idx, self.mapobj.res_arr, smax_lf
                    )
                linefit.step = self.step
                alpha = linefit.calc_fval_for_different_kvalues_at_this_step(self.e1)
                end_lf = timer()
                if timeit:
                    print(" time for line fit: ", end_lf - start_lf)
                # new linefit
                """ if i == 0:
                    lft = linefit_class.linefit2()
                lft.cbin_idx = self.mapobj.cbin_idx
                lft.cbin = self.mapobj.cbin
                lft.res_arr = self.mapobj.res_arr
                lft.smax = smax_lf
                lft.e0 = self.e0
                lft.e1 = self.e1
                lft.get_linefit_static_data()
                lft.step = self.step
                #alpha_t = lft.scalar_opt_trans()
                lft.q_prev = self.q
                #alpha_r = lft.scalar_opt_rot()
                lft.scalar_opt_rot()
                alpha_t = lft.alpha_t
                alpha_r = lft.alpha_r
                self.t = self.step[:3] * alpha_t
                tmp = np.insert(self.step[3:] * alpha_r, 0, 0.0) """
                # translation
                self.t = self.step[:3] * alpha[0]
                t_accum = t_accum + self.t
                translation_vec = trans_in_angstrom(t_accum, self.pixsize, self.ful_dim)
                # rotation
                tmp = np.insert(self.step[3:] * alpha[1], 0, 0.0)
                q_accum = q_accum + tmp
                q_accum = q_accum / np.sqrt(np.dot(q_accum, q_accum))
                rm_accum = quaternions.get_RM(q_accum)
                theta2 = np.arccos((np.trace(rm_accum) - 1) / 2) * 180.0 / np.pi
                tmp = tmp + q_init
                q = tmp / np.sqrt(np.dot(tmp, tmp))
                self.q = q
                fval_previous = fval
                end = timer()
                if timeit:
                    print("time for one cycle:", end - start)


def trans_in_angstrom(t, pixsize, dim):
    t_angs = np.asarray(t) * np.asarray(pixsize) * np.asarray(dim)
    return np.sqrt(np.sum(t_angs * t_angs))
