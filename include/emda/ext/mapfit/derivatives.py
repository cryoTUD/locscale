"""
Author: "Rangana Warshamanage, Garib N. Murshudov"
MRC Laboratory of Molecular Biology
    
This software is released under the
Mozilla Public License, version 2.0; see LICENSE.
"""

from __future__ import absolute_import, division, print_function, unicode_literals
import numpy as np
import math
import sys
import fcodes_fast
from emda.core.quaternions import derivatives_wrt_q, quart2axis
from timeit import default_timer as timer

timeit = False


def new_dFs2(mapin, xyz, vol):
    nx, ny, nz = mapin.shape
    dfs = np.zeros(shape=(nx, ny, nz, 3), dtype=np.complex64)
    for i in range(3):
        dfs[:, :, :, i] = np.fft.fftshift(
            (1 / vol) * 2j * np.pi * np.fft.fftn(mapin * xyz[i])
        )
    return dfs


def get_dqda(q):
    _,_,_,angle = quart2axis(q)
    dqda = np.zeros(3, dtype='float')
    for i in range(3):
        dqda[i] = math.sin(angle/2)
    return dqda


def new_derivatives(e0, e1, w_grid, w2_grid, q, sv, xyz, xyz_sum, vol, dfrs=None):
    nx, ny, nz = e0.shape
    start = timer()
    sv_np = np.zeros((nx, ny, nz, 3), dtype=np.float64)
    for i in range(3):
        sv_np[:, :, :, i] = sv[i]
    dRdq = derivatives_wrt_q(q)
    if dfrs is None:
        start = timer()
        dFRs = new_dFs2(np.real(np.fft.ifftn(np.fft.ifftshift(e1))), xyz, vol)
        end = timer()
        if timeit:
            print(" time for dFRs calculation: ", end - start)
    else:
        print("Interpolated DFRs are used!")
        dFRs = dfrs
    start = timer()
    df_val, ddf_val = fcodes_fast.calc_derivatives(
        e0, e1, w_grid, w2_grid, sv_np, dFRs, dRdq, xyz_sum, vol, nx, ny, nz
    )
    end = timer()
    if timeit:
        print(" time for derivative calculation: ", end - start)
    if np.linalg.cond(ddf_val) < 1 / sys.float_info.epsilon:
        ddf_val_inv = np.linalg.pinv(ddf_val)
        step = ddf_val_inv.dot(-df_val)
        return step, df_val
    else:
        raise SystemExit("Derivative matrix is non-invertible! Stopping now...")


def new_derivatives2(e0, e1, w_grid, w2_grid, q, sv, xyz, xyz_sum, vol, dfrs=None):
    nx, ny, nz = e0.shape
    start = timer()
    sv_np = np.zeros((nx, ny, nz, 3), dtype=np.float64)
    for i in range(3):
        sv_np[:, :, :, i] = sv[i]
    dRdq = derivatives_wrt_q(q)
    dqda = get_dqda(q)
    dqda_x = np.array([1.0, 1.0, 1.0], dtype='float')
    if dfrs is None:
        start = timer()
        dFRs = new_dFs2(np.real(np.fft.ifftn(np.fft.ifftshift(e1))), xyz, vol)
        end = timer()
        if timeit:
            print(" time for dFRs calculation: ", end - start)
    else:
        print("Interpolated DFRs are used!")
        dFRs = dfrs
    start = timer()
    df_val, ddf_val = fcodes_fast.calc_derivatives2(
        e0, e1, w_grid, w2_grid, sv_np, dFRs, dRdq, dqda_x, xyz_sum, vol, nx, ny, nz
    )
    end = timer()
    if timeit:
        print(" time for derivative calculation: ", end - start)
    if np.linalg.cond(ddf_val) < 1 / sys.float_info.epsilon:
        ddf_val_inv = np.linalg.pinv(ddf_val)
        step = ddf_val_inv.dot(-df_val)
        # calculating step on axes
        df_ax = np.zeros(3, dtype='float')
        ddf_ax = np.zeros((3,3), dtype='float')
        for i in range(len(dqda)):
            df_ax[i] = df_val[3+i] * dqda[i]
            for j in range(len(dqda)):
                ddf_ax[i,j] = ddf_val[3+i,3+j] * dqda[i] * dqda[j]
        ddf_ax_inv = np.linalg.pinv(ddf_ax)
        step_ax = ddf_ax_inv.dot(-df_ax)
        return step, step_ax
    else:
        raise SystemExit("Derivative matrix is non-invertible! Stopping now...")


def axis_derivatives(e0, e1, w_grid, w2_grid, q, sv, xyz, xyz_sum, vol, dfrs=None, axis=None):
    nx, ny, nz = e0.shape
    start = timer()
    sv_np = np.zeros((nx, ny, nz, 3), dtype=np.float64)
    for i in range(3):
        sv_np[:, :, :, i] = sv[i]
    dRdq = derivatives_wrt_q(q)
    az, ay, ax = axis[0], axis[1], axis[2]
    _, _, _, angle = quart2axis(q)
    dqda = get_dqda(q)
    if dfrs is None:
        # print('Calculating derivatives...')
        start = timer()
        dFRs = new_dFs2(np.real(np.fft.ifftn(np.fft.ifftshift(e1))), xyz, vol)
        end = timer()
        if timeit:
            print(" time for dFRs calculation: ", end - start)
    else:
        print("Interpolated DFRs are used!")
        dFRs = dfrs
    start = timer()
    df_val, ddf_val = fcodes_fast.calc_derivatives(
        e0, e1, w_grid, w2_grid, sv_np, dFRs, dRdq, xyz_sum, vol, nx, ny, nz
    )
    end = timer()
    if timeit:
        print(" time for derivative calculation: ", end - start)
    if np.linalg.cond(ddf_val) < 1 / sys.float_info.epsilon:
        ddf_val_inv = np.linalg.pinv(ddf_val)
    else:
        raise SystemExit("Derivative matrix is non-invertible! Stopping now...")
    df_daz = df_val[3] * dqda[0] + df_val[5] * (az/ax) * dqda[2]
    df_day = df_val[4] * dqda[1] + df_val[5] * (ay/ax) * dqda[2]
    ddf_zz = (ddf_val[3,3] - 2 * ddf_val[3,5] * (az/ax) + ddf_val[5,5] * (az/ax)*(az/ax)) * dqda[0] * dqda[0]
    ddf_yy = (ddf_val[4,4] - 2 * ddf_val[4,5] * (ay/ax) + ddf_val[5,5] * (ay/ax)*(ay/ax)) * dqda[1] * dqda[1]
    ddf_zy = ddf_val[3,4] * (az/ax) * (ay/ax) * dqda[0] * dqda[1]
    ddf_yz = ddf_val[4,3] * (az/ax) * (ay/ax) * dqda[1] * dqda[0]

    df_axis = np.zeros(2, dtype='float')
    ddf_axis = np.zeros((2,2), dtype='float')
    df_axis[0] = df_daz
    df_axis[1] = df_day
    ddf_axis[0,0] = ddf_zz
    ddf_axis[0,1] = ddf_zy
    ddf_axis[1,0] = ddf_yz
    ddf_axis[1,1] = ddf_yy
    ddf_axis_inv = np.linalg.pinv(ddf_axis)
    step_axis = ddf_axis_inv.dot(-df_axis)
    step = ddf_val_inv.dot(-df_val)
    return step, step_axis

