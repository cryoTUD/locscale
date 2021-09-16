"""
Author: "Rangana Warshamanage, Garib N. Murshudov"
MRC Laboratory of Molecular Biology
    
This software is released under the
Mozilla Public License, version 2.0; see LICENSE.
"""
from __future__ import absolute_import, division, print_function, unicode_literals
import numpy as np
from emda.core import restools, quaternions, fsc
import fcodes_fast
from emda.ext.mapfit.utils import get_FRS

np.set_printoptions(suppress=True)  # Suppress insignificant values for clarity


class LineFit:
    def __init__(self):
        self.e0_lf = None
        self.e1_lf = None
        self.cbin_idx_lf = None
        self.cbin_lf = None
        self.w_grid = None
        self.step = []

    def get_linefit_static_data(self, e0, cbin_idx, res_arr, smax):
        (
            self.e0_lf,
            self.cbin_idx_lf,
            self.cbin_lf,
        ) = restools.cut_resolution_for_linefit(e0, cbin_idx, res_arr, smax)

    def f(self, k):
        from emda.ext.mapfit.utils import get_fsc_wght, get_FRS

        # w = 1.0 for line search
        nx, ny, nz = self.e0_lf.shape
        t = np.asarray(self.step[:3], 'float') * k[0]
        st, _, _, _ = fcodes_fast.get_st(nx, ny, nz, t)
        q_init = np.array([1.0, 0.0, 0.0, 0.0])
        tmp = np.insert(np.asarray(self.step[3:], 'float') * k[1], 0, 0.0)
        tmp = tmp + q_init
        q = tmp / np.sqrt(np.dot(tmp, tmp))
        rotmat = quaternions.get_RM(q)
        #ers = get_FRS(rotmat, self.e1_lf * st, interp="linear")
        ers = get_FRS(rotmat, self.e1_lf, interp="linear")
        w_grid = get_fsc_wght(
            self.e0_lf, ers[:, :, :, 0] * st, self.cbin_idx_lf, self.cbin_lf
        )
        fval = np.sum(w_grid * self.e0_lf * np.conjugate(ers[:, :, :, 0] * st))
        return -fval.real

    def calc_fval_for_different_kvalues_at_this_step(self, e1):
        from scipy import optimize

        nx, ny, nz = e1.shape
        cx = self.e0_lf.shape[0] // 2
        cy = cz = cx
        dx = int((nx - 2 * cx) / 2)
        dy = int((ny - 2 * cy) / 2)
        dz = int((nz - 2 * cz) / 2)
        self.e1_lf = e1[dx : dx + 2 * cx, dy : dy + 2 * cy, dz : dz + 2 * cz]
        assert self.e1_lf.shape == self.e0_lf.shape
        init_guess = np.array([1.0, 1.0], dtype='float')
        minimum = optimize.minimize(self.f, init_guess, method="Powell")
        return minimum.x


class linefit2:
    def __init__(self):
        self.e0 = None
        self.e1 = None
        self.cbin_idx = None
        self.cbin = None
        self.res_arr = None
        self.smax = None
        self.step = None
        self.t_ini = None
        self.q_prev = None
        self.alpha_t = None
        self.alpha_r = None

    def get_linefit_static_data(self):
        (
            self.e0,
            self.bin_idx,
            self.nbin,
        ) = restools.cut_resolution_for_linefit(self.e0, self.cbin_idx, self.res_arr, self.smax)
        (    
            self.e1,
            _,
            _,
        ) = restools.cut_resolution_for_linefit(self.e1, self.cbin_idx, self.res_arr, self.smax)


    def get_fsc_wght(self, e0, ert, bin_idx, nbin):
        cx, cy, cz = e0.shape
        bin_stats = fsc.anytwomaps_fsc_covariance(e0, ert, bin_idx, nbin)
        bin_fsc, _ = bin_stats[0], bin_stats[1]
        bin_fsc = np.array(bin_fsc, dtype=np.float64, copy=False)
        w_grid = fcodes_fast.read_into_grid(bin_idx, bin_fsc, nbin, cx, cy, cz)
        w_grid = np.array(w_grid, dtype=np.float64, copy=False)
        return w_grid

    def func_t(self, i):
        nx, ny, nz = self.e0.shape
        t = self.step[:3] * i
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

    def func_r(self, i):
        tmp = np.insert(self.step[3:] * i, 0, 0.0)
        q = tmp + self.q_prev
        q = q / np.sqrt(np.dot(q, q))
        rotmat = quaternions.get_RM(q)
        ers = get_FRS(rotmat, self.e1, interp="linear")
        w_grid = self.get_fsc_wght(self.e0, ers[:, :, :, 0], self.bin_idx, self.nbin)
        fval = np.real(np.sum(w_grid * self.e0 * np.conjugate(ers[:, :, :, 0])))
        return -fval

    def scalar_opt_rot(self, t=None):
        from scipy.optimize import minimize_scalar

        nx, ny, nz = self.e1.shape
        """ self.alpha_t = self.scalar_opt_trans()
        t = self.alpha_t * self.step[:3]
        if t is not None:
            st, _, _, _ = fcodes_fast.get_st(nx, ny, nz, t)
            self.e1 = self.e1 * st """
        f = self.func_r
        res = minimize_scalar(f, method="brent")
        self.alpha_r = res.x
        return res.x