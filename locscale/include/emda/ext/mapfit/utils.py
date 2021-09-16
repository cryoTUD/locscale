"""
Author: "Rangana Warshamanage, Garib N. Murshudov"
MRC Laboratory of Molecular Biology
    
This software is released under the
Mozilla Public License, version 2.0; see LICENSE.
"""

from __future__ import absolute_import, division, print_function, unicode_literals

import numpy as np
import fcodes_fast
from emda import core, ext
from timeit import default_timer as timer
from emda.config import debug_mode


def get_FRS(RM, E2, interp='linear'):
    if len(E2.shape) == 3:
        E2 = np.expand_dims(E2, axis=3)
    ERS = get_interp(RM, E2, interp)
    return ERS


def get_interp(RM, data, interp='linear'):
    assert len(data.shape) == 4
    ih, ik, il, n = data.shape
    if interp == "cubic":
        interp3d = fcodes_fast.tricubic(RM, data, debug_mode, n, ih, ik, il)
    if interp == "linear":
        interp3d = fcodes_fast.trilinear(RM, data, debug_mode, n, ih, ik, il)
    return interp3d


def create_xyz_grid(uc, nxyz):
    x = np.fft.fftfreq(nxyz[0]) #* uc[0]
    y = np.fft.fftfreq(nxyz[1]) #* uc[1]
    z = np.fft.fftfreq(nxyz[2]) #* uc[2]
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


def avg_and_diffmaps(
    maps2avg,
    uc,
    nbin,
    sgnl_var,
    totl_var,
    covar,
    hffsc,
    bin_idx,
    s_grid,
    res_arr,
    Bf_arr,
):
    # average and difference map calculation
    import numpy.ma as ma

    nx, ny, nz = maps2avg[0].shape
    nmaps = len(maps2avg)
    unit_cell = uc
    all_maps = np.zeros(shape=(nx, ny, nz, nmaps), dtype="complex")
    for i in range(nmaps):
        all_maps[:, :, :, i] = maps2avg[i]
    print(all_maps.shape)
    #
    S_mat = np.zeros(shape=(nmaps, nmaps, nbin), dtype="float")
    T_mat = np.zeros(shape=(nmaps, nmaps, nbin), dtype="float")
    F_mat = np.zeros(shape=(nmaps, nmaps, nbin), dtype="float")
    # Populate Sigma matrix
    start = timer()
    for i in range(nmaps):
        for j in range(nmaps):
            if i == j:
                print("i,j=", i, j)
                # Diagonals
                S_mat[i, i, :] = 1.0
                T_mat[i, i, :] = 1.0
                F_mat[i, i, :] = np.sqrt(
                    ma.masked_less_equal(hffsc[i], 0.0).filled(0.0)
                )
            elif i == 0 and j > i:
                # Off diagonals
                print("i,j=", i, j)
                core.plotter.plot_nlines_log(
                    res_arr,
                    [covar[j - 1], sgnl_var[i], sgnl_var[j]],
                    ["S12", "S11", "S22"],
                    "log_variance_signal.eps",
                )
                Sigma_ij_ma = ma.masked_less_equal(covar[j - 1], 0.0)
                sgnl_var_ij_ma = ma.masked_where(
                    ma.getmask(Sigma_ij_ma), sgnl_var[i] * sgnl_var[j]
                )
                S_mat[i, j, :] = (Sigma_ij_ma / np.sqrt(sgnl_var_ij_ma)).filled(0.0)
                Totalv_ij_ma = ma.masked_where(
                    ma.getmask(Sigma_ij_ma), totl_var[i] * totl_var[j]
                )
                T_mat[i, j, :] = (Sigma_ij_ma / np.sqrt(Totalv_ij_ma)).filled(0.0)
                core.plotter.plot_nlines_log(
                    res_arr,
                    [covar[j - 1], totl_var[i], totl_var[j]],
                    ["S12", "T11", "T22"],
                    "log_variance_totalv.eps",
                )
            elif j == 0 and i > j:
                print("i,j=", i, j)
                S_mat[i, j, :] = S_mat[j, i, :]
                T_mat[i, j, :] = T_mat[j, i, :]
            # else:
            #    print('i,j=',i,j)
            #    S_mat[i,j,:] = S_mat[j,i,:]
            #    T_mat[i,j,:] = T_mat[j,i,:]
    # Plotting
    core.plotter.plot_nlines(
        res_arr,
        [S_mat[0, 0, :], S_mat[0, 1, :], S_mat[1, 0, :], S_mat[1, 1, :]],
        "S_mat_fsc_ij.eps",
        ["FSC11", "FSC12", "FSC21", "FSC22"],
    )
    core.plotter.plot_nlines(
        res_arr,
        [F_mat[0, 0, :], F_mat[1, 1, :]],
        "F_mat_fsc_ij.eps",
        ["sqrt(FSC11)", "sqrt(FSC22)"],
    )
    core.plotter.plot_nlines(
        res_arr,
        [T_mat[0, 0, :], T_mat[0, 1, :], T_mat[1, 0, :], T_mat[1, 1, :]],
        "T_mat_fsc_ij.eps",
        ["FSC11", "FSC12", "FSC21", "FSC22"],
    )
    end = timer()
    print("Time for sigma matrix population: ", end - start)
    start = timer()
    # Variance weighted matrices calculation
    Wgt = np.zeros(shape=(nmaps, nmaps, nbin))
    for ibin in range(nbin):
        T_mat_inv = np.linalg.pinv(T_mat[:, :, ibin])  # Moore-Penrose psedo-inversion
        tmp = np.dot(F_mat[:, :, ibin], T_mat_inv)
        Wgt[:, :, ibin] = np.dot(S_mat[:, :, ibin], tmp)
    core.plotter.plot_nlines(
        res_arr,
        [Wgt[0, 0, :], Wgt[0, 1, :], Wgt[1, 0, :], Wgt[1, 1, :]],
        "Wgt_map_ij.eps",
        ["W11", "W12", "W21", "W22"],
    )
    end = timer()
    print("time for pinv  ", end - start)
    # output data
    fsmat = open("smat.txt", "w")
    ftmat = open("tmat.txt", "w")
    # ftmatinv = open("tmatinv.txt", "w")
    fwmat = open("wmat.txt", "w")
    for i in range(nbin):
        fsmat.write(
            "{:.2f} {:.4f} {:.4f} {:.4f} {:.4f}\n".format(
                res_arr[i],
                S_mat[0, 0, i],
                S_mat[0, 1, i],
                S_mat[1, 0, i],
                S_mat[1, 1, i],
            )
        )
        ftmat.write(
            "{:.2f} {:.4f} {:.4f} {:.4f} {:.4f}\n".format(
                res_arr[i],
                T_mat[0, 0, i],
                T_mat[0, 1, i],
                T_mat[1, 0, i],
                T_mat[1, 1, i],
            )
        )
        fwmat.write(
            "{:.2f} {:.4f} {:.4f} {:.4f} {:.4f}\n".format(
                res_arr[i], Wgt[0, 0, i], Wgt[0, 1, i], Wgt[1, 0, i], Wgt[1, 1, i]
            )
        )

    # Average map calculation
    nbf = len(Bf_arr)
    assert all_maps.shape[:3] == bin_idx.shape
    assert all_maps.shape[3] == nmaps
    AVG_Maps = fcodes_fast.calc_avg_maps(
        all_maps,
        bin_idx,
        s_grid,
        Wgt,
        Bf_arr,
        unit_cell,
        debug_mode,
        nbin,
        nmaps,
        nbf,
        nx,
        ny,
        nz,
    )
    return AVG_Maps


def output_maps(
    averagemaps, com_lst, t_lst, r_lst, center, unit_cell, map_origin, bf_arr
):

    start = timer()
    nx, ny, nz = averagemaps.shape[:3]
    # center = (nx/2, ny/2, nz/2)
    for ibf in range(averagemaps.shape[4]):
        if bf_arr[ibf] < 0.0:
            Bcode = "_blur" + str(abs(bf_arr[ibf]))
        elif bf_arr[ibf] > 0.0:
            Bcode = "_sharp" + str(abs(bf_arr[ibf]))
        else:
            Bcode = "_unsharpened"
        for imap in range(averagemaps.shape[3]):
            filename_mrc = "avgmap_" + str(imap) + Bcode + ".mrc"
            if com_lst:
                com = com_lst[imap]
            if imap > 0:
                nx, ny, nz = averagemaps[:, :, :, imap, ibf].shape
                if r_lst and t_lst:
                    rotmat = r_lst[imap - 1]
                    avgmap = get_FRS(
                        np.transpose(rotmat),
                        averagemaps[:, :, :, imap, ibf],
                        interp="cubic",
                    )[:, :, :, 0]
                    t = -1.0 * np.asarray(t_lst[imap - 1], dtype="float")  # inverse t
                    st, _, _, _ = fcodes_fast.get_st(nx, ny, nz, t)
                    if com_lst and (center is not None):
                        data2write = apply_shift(
                            f_map=avgmap * st, com=com, center=center
                        )
                    else:
                        data2write = apply_shift(f_map=avgmap)
                else:
                    data2write = apply_shift(f_map=averagemaps[:, :, :, imap, ibf])
            if imap == 0:
                if com_lst and (center is not None):
                    data2write = apply_shift(
                        f_map=averagemaps[:, :, :, imap, ibf], com=com, center=center
                    )
                else:
                    data2write = apply_shift(f_map=averagemaps[:, :, :, imap, ibf])
            core.iotools.write_mrc(data2write, filename_mrc, unit_cell, map_origin)

        ## Difference and Biase removed map calculation
        # nmaps = averagemaps.shape[3]
        # for m in range(nmaps-1):
        #    for n in range(m+1,nmaps):
        #        fname_diff1 = 'diffmap_m'+str(n)+'-m'+str(m)+Bcode
        #        com = com_lst[n]
        #        dm1 = apply_shift(averagemaps[:,:,:,n,ibf] - averagemaps[:,:,:,m,ibf],com,center)
        #        iotools.write_mrc(dm1,
        #                          fname_diff1+'.mrc',
        #                          unit_cell,map_origin)
        #        '''write_3d2mtz_fortran(unit_cell,dm1,fname_diff1+'.mtz')'''
        #        dm2 = apply_shift(averagemaps[:,:,:,m,ibf] - averagemaps[:,:,:,n,ibf],com,center)
        #        fname_diff2 = 'diffmap_m'+str(m)+'-m'+str(n)+Bcode
        #        iotools.write_mrc(dm2,
        #                          fname_diff2+'.mrc',
        #                          unit_cell,map_origin)
        #        '''write_3d2mtz_fortran(unit_cell,dm2,fname_diff2+'.mtz')'''
        ## 2Fo-Fc type maps
        # twom2m1 = apply_shift(averagemaps[:,:,:,n,ibf] + dm1, com, center)
        # fname1_2FoFc = '2m'+str(n)+'-m'+str(m)
        # iotools.write_mrc(twom2m1,
        #                  fname1_2FoFc+'.mrc',
        #                  unit_cell,map_origin)
        #'''write_3d2mtz_fortran(unit_cell,twom2m1,fname1_2FoFc+'.mtz')'''
        # twom1m2 = apply_shift(averagemaps[:,:,:,m,ibf] + dm2, com, center)
        # fname2_2FoFc = '2m'+str(m)+'-m'+str(n)
        # iotools.write_mrc(twom1m2,
        #                  fname2_2FoFc+'.mrc',
        #                  unit_cell,map_origin)
        #'''write_3d2mtz_fortran(unit_cell,twom1m2,fname2_2FoFc+'.mtz')'''
    end = timer()
    print("Map output time: ", end - start)


def apply_shift(f_map, com=None, center=None):
    from scipy.ndimage.interpolation import shift

    data = np.real(np.fft.ifftshift(np.fft.ifftn(np.fft.ifftshift(f_map))))
    if com is not None:
        if center is not None:
            data = shift(data, np.subtract(com, center))
    return data


def remove_unwanted_corners(uc, target_dim):
    nx, ny, nz = target_dim
    fResArr = core.restools.get_resArr(uc, nx)
    cut_mask = core.restools.remove_edge(fResArr, fResArr[-1])
    cc_mask = np.zeros(shape=(nx, ny, nz), dtype="int")
    cx, cy, cz = cut_mask.shape
    dx = (nx - cx) // 2
    dy = (ny - cy) // 2
    dz = (nz - cz) // 2
    print(dx, dy, dz)
    cc_mask[dx : dx + cx, dy : dy + cy, dz : dz + cz] = cut_mask
    return cc_mask


def double_the_axes(arr1):
    nx, ny, nz = arr1.shape
    big_arr1 = np.zeros((2 * nx, 2 * ny, 2 * nz), dtype=arr1.dtype)
    dx = int(nx / 2)
    dy = int(ny / 2)
    dz = int(nz / 2)
    big_arr1[dx : dx + nx, dy : dy + ny, dz : dz + nz] = arr1
    return big_arr1


def get_minimum_bounding_dims(arr1):
    from scipy import ndimage

    # box dim check begin
    indices = np.indices(arr1.shape)
    map_sum = np.sum(arr1) / 1000
    arr1 = arr1 * (arr1 > map_sum)
    com = ndimage.measurements.center_of_mass(arr1)
    dist2 = (
        (indices[0, :, :] - com[0]) ** 2
        + (indices[1, :, :] - com[1]) ** 2
        + (indices[2, :, :] - com[2]) ** 2
    )
    # np.max(dist2[mask > 0.5])
    min_dim = int(np.sqrt(np.max(dist2)) * 2) + 1
    print(com, min_dim)
    # box dim check end
    return com, min_dim


def realsp_map_interpolation(arr1, RM):
    from scipy import ndimage

    nx, ny, nz = arr1.shape
    rotated_map = fcodes_fast.trilinear_map(RM, arr1, debug_mode, nx, ny, nz)
    return rotated_map


def calc_averagemaps_simple(
    maps2avg,
    uc,
    nbin,
    sgnl_var,
    totl_var,
    covar,
    bin_idx,
    s_grid,
    res_arr,
    Bf_arr,
    com_lst,
    box_centr,
    origin,
):
    # average and difference map calculation

    nx, ny, nz = maps2avg[0].shape
    nmaps = len(maps2avg)
    unit_cell = uc
    all_maps = np.zeros(shape=(nx, ny, nz, nmaps), dtype="complex")
    for i in range(nmaps):
        all_maps[:, :, :, i] = maps2avg[i]
    print(all_maps.shape)
    #
    S_mat = np.zeros(shape=(nmaps, nmaps, nbin), dtype="float")  # signal variance
    T_mat = np.zeros(shape=(nmaps, nmaps, nbin), dtype="float")  # total variance
    # Populate Sigma matrix
    start = timer()
    for i in range(nmaps):
        for j in range(nmaps):
            if i == j:
                print("i,j=", i, j)
                # Diagonals
                S_mat[i, i, :] = sgnl_var[i]
                T_mat[i, i, :] = totl_var[i]
            elif i == 0 and j > i:
                # Off diagonals
                print("i,j=", i, j)
                S_mat[i, j, :] = covar[j - 1]
                T_mat[i, j, :] = covar[j - 1]
            elif j == 0 and i > j:
                print("i,j=", i, j)
                S_mat[i, j, :] = S_mat[j, i, :]
                T_mat[i, j, :] = T_mat[j, i, :]
            # else:
            #    print('i,j=',i,j)
            #    S_mat[i,j,:] = S_mat[j,i,:]
            #    T_mat[i,j,:] = T_mat[j,i,:]
    # Plotting
    core.plotter.plot_nlines_log(
        res_arr,
        [S_mat[0, 0, :], S_mat[0, 1, :], S_mat[1, 0, :], S_mat[1, 1, :]],
        ["S11", "S12", "S21", "S22"],
        "S_mat_fsc_ij.eps",
    )
    core.plotter.plot_nlines_log(
        res_arr,
        [T_mat[0, 0, :], T_mat[0, 1, :], T_mat[1, 0, :], T_mat[1, 1, :]],
        ["T11", "T12", "T21", "T22"],
        "T_mat_fsc_ij.eps",
    )

    end = timer()
    print("Time for sigma matrix population: ", end - start)
    start = timer()
    # Variance weighted matrices calculation
    Wgt = np.zeros(shape=(nmaps, nmaps, nbin))
    for ibin in range(nbin):
        T_mat_inv = np.linalg.pinv(T_mat[:, :, ibin])  # Moore-Penrose psedo-inversion
        Wgt[:, :, ibin] = np.dot(S_mat[:, :, ibin], T_mat_inv)
    core.plotter.plot_nlines(
        res_arr,
        [Wgt[0, 0, :], Wgt[0, 1, :], Wgt[1, 0, :], Wgt[1, 1, :]],
        "Wgt_map_ij.eps",
        ["W11", "W12", "W21", "W22"],
    )
    end = timer()
    print("time for pinv  ", end - start)
    # output data
    fsmat = open("smat.txt", "w")
    ftmat = open("tmat.txt", "w")
    # ftmatinv = open("tmatinv.txt", "w")
    fwmat = open("wmat.txt", "w")
    for i in range(nbin):
        fsmat.write(
            "{:.2f} {:.4f} {:.4f} {:.4f} {:.4f}\n".format(
                res_arr[i],
                S_mat[0, 0, i],
                S_mat[0, 1, i],
                S_mat[1, 0, i],
                S_mat[1, 1, i],
            )
        )
        ftmat.write(
            "{:.2f} {:.4f} {:.4f} {:.4f} {:.4f}\n".format(
                res_arr[i],
                T_mat[0, 0, i],
                T_mat[0, 1, i],
                T_mat[1, 0, i],
                T_mat[1, 1, i],
            )
        )
        fwmat.write(
            "{:.2f} {:.4f} {:.4f} {:.4f} {:.4f}\n".format(
                res_arr[i], Wgt[0, 0, i], Wgt[0, 1, i], Wgt[1, 0, i], Wgt[1, 1, i]
            )
        )

    # Average map calculation
    nbf = len(Bf_arr)
    assert all_maps.shape[:3] == bin_idx.shape
    assert all_maps.shape[3] == nmaps
    from scipy.ndimage.interpolation import shift

    """for i in range(all_maps.shape[3]):
        data2write = np.real(np.fft.ifftshift(np.fft.ifftn(np.fft.ifftshift(all_maps[:,:,:,i]))))
        data2write = shift(data2write, np.subtract(com1,box_centr))
        fname = 'allmaps_'+str(i)+'.mrc'
        iotools.write_mrc(data2write,fname,unit_cell,origin)"""
    AVG_Maps = fcodes_fast.calc_avg_maps(
        all_maps,
        bin_idx,
        s_grid,
        Wgt,
        Bf_arr,
        unit_cell,
        debug_mode,
        nbin,
        nmaps,
        nbf,
        nx,
        ny,
        nz,
    )
    print("avg_map shape: ", AVG_Maps.shape)
    """for i in range(AVG_Maps.shape[3]):
        data2write = np.real(np.fft.ifftshift(np.fft.ifftn(np.fft.ifftshift(AVG_Maps[:,:,:,i,0]))))
        fname1 = 'avgmaps_before'+str(i)+'.mrc'
        iotools.write_mrc(data2write,fname1,unit_cell,origin)
        data2write = shift(data2write, np.subtract(com1,box_centr))
        fname = 'avgmaps_'+str(i)+'.mrc'
        iotools.write_mrc(data2write,fname,unit_cell,origin)"""
    return AVG_Maps


def calc_averagemaps_simple_allmaps(
    maps2avg, uc, nbin, varcovmat, totl_var, bin_idx, s_grid, res_arr, Bf_arr
):
    # average and difference map calculation

    nx, ny, nz = maps2avg[0].shape
    nmaps = len(maps2avg)
    unit_cell = uc
    all_maps = np.zeros(shape=(nx, ny, nz, nmaps), dtype="complex")
    for i in range(nmaps):
        all_maps[:, :, :, i] = maps2avg[i]
    print(all_maps.shape)
    #
    S_mat = np.zeros(shape=(nmaps, nmaps, nbin), dtype="float")  # signal variance
    T_mat = np.zeros(shape=(nmaps, nmaps, nbin), dtype="float")  # total variance
    # Populate Sigma matrix
    start = timer()
    for i in range(nmaps):
        for j in range(nmaps):
            if i == j:
                print("i,j=", i, j)
                # Diagonals
                S_mat[i, i, :] = varcovmat[i, i, :]
                T_mat[i, i, :] = totl_var[i]
            elif j > i:
                # Off diagonals
                print("i,j=", i, j)
                S_mat[i, j, :] = varcovmat[i, j, :]
                T_mat[i, j, :] = varcovmat[i, j, :]
            elif i > j:
                print("i,j=", i, j)
                S_mat[i, j, :] = varcovmat[j, i, :]
                T_mat[i, j, :] = varcovmat[j, i, :]
    # Plotting
    core.plotter.plot_nlines_log(
        res_arr,
        [S_mat[0, 0, :], S_mat[0, 1, :], S_mat[1, 0, :], S_mat[1, 1, :]],
        ["S11", "S12", "S21", "S22"],
        "S_mat_fsc_ij.eps",
    )
    core.plotter.plot_nlines_log(
        res_arr,
        [T_mat[0, 0, :], T_mat[0, 1, :], T_mat[1, 0, :], T_mat[1, 1, :]],
        ["T11", "T12", "T21", "T22"],
        "T_mat_fsc_ij.eps",
    )

    end = timer()
    print("Time for sigma matrix population: ", end - start)
    start = timer()
    # Variance weighted matrices calculation
    ftmatinv = open("tmatinv.txt", "w")
    Wgt = np.zeros(shape=(nmaps, nmaps, nbin))
    for ibin in range(nbin):
        T_mat_inv = np.linalg.pinv(T_mat[:, :, ibin])  # Moore-Penrose psedo-inversion
        ftmatinv.write(
            "{:.2f} {:.4f} {:.4f} {:.4f} {:.4f}\n".format(
                res_arr[ibin],
                T_mat_inv[0, 0],
                T_mat_inv[0, 1],
                T_mat_inv[1, 0],
                T_mat_inv[1, 1],
            )
        )
        # Wgt[:,:,ibin] = np.dot(S_mat[:,:,ibin],T_mat_inv)
        Wgt[:, :, ibin] = np.matmul(S_mat[:, :, ibin], T_mat_inv)
    core.plotter.plot_nlines(
        res_arr,
        [Wgt[0, 0, :], Wgt[0, 1, :], Wgt[1, 0, :], Wgt[1, 1, :]],
        "Wgt_map_ij.eps",
        ["W11", "W12", "W21", "W22"],
    )
    end = timer()
    print("time for pinv  ", end - start)
    # output data
    fsmat = open("smat.txt", "w")
    ftmat = open("tmat.txt", "w")
    # ftmatinv = open("tmatinv.txt", "w")
    fwmat = open("wmat.txt", "w")
    for i in range(nbin):
        fsmat.write(
            "{:.2f} {:.4f} {:.4f} {:.4f} {:.4f}\n".format(
                res_arr[i],
                S_mat[0, 0, i],
                S_mat[0, 1, i],
                S_mat[1, 0, i],
                S_mat[1, 1, i],
            )
        )
        ftmat.write(
            "{:.2f} {:.4f} {:.4f} {:.4f} {:.4f}\n".format(
                res_arr[i],
                T_mat[0, 0, i],
                T_mat[0, 1, i],
                T_mat[1, 0, i],
                T_mat[1, 1, i],
            )
        )
        fwmat.write(
            "{:.2f} {:.4f} {:.4f} {:.4f} {:.4f}\n".format(
                res_arr[i], Wgt[0, 0, i], Wgt[0, 1, i], Wgt[1, 0, i], Wgt[1, 1, i]
            )
        )

    # Average map calculation
    nbf = len(Bf_arr)
    assert all_maps.shape[:3] == bin_idx.shape
    assert all_maps.shape[3] == nmaps
    AVG_Maps = fcodes_fast.calc_avg_maps(
        all_maps,
        bin_idx,
        s_grid,
        Wgt,
        Bf_arr,
        unit_cell,
        debug_mode,
        nbin,
        nmaps,
        nbf,
        nx,
        ny,
        nz,
    )
    print("avg_map shape: ", AVG_Maps.shape)
    return AVG_Maps


def set_dim_even(x):
    # Check if dims are even
    """ if x.shape[0] % 2 != 0:
        xshape = list(x.shape)
        xshape[0] = xshape[0] + 1
        xshape[1] = xshape[1] + 1
        xshape[2] = xshape[2] + 1
        temp = np.zeros(xshape, x.dtype)
        temp[:-1, :-1, :-1] = x
        x = temp """
    xshape = list(x.shape)
    xshape = (np.array(xshape) + np.array(xshape) % 2).tolist()
    temp = np.zeros(xshape, x.dtype)
    temp[-x.shape[0]:, -x.shape[1]:, -x.shape[2]:] = x
    x = temp
    return x


def get_fsc_wght(e0, ert, bin_idx, nbin):
    cx, cy, cz = e0.shape
    bin_stats = core.fsc.anytwomaps_fsc_covariance(e0, ert, bin_idx, nbin)
    fsc, f1f2_covar = bin_stats[0], bin_stats[1]
    fsc = np.array(fsc, dtype=np.float64, copy=False)
    w_grid = fcodes_fast.read_into_grid(bin_idx, fsc / (1 - fsc ** 2), nbin, cx, cy, cz)
    w_grid = np.array(w_grid, dtype=np.float64, copy=False)
    return w_grid


def frac2ang(t_frac, pix_size):
    t_ang = np.asarray(t_frac, dtype="float") * pix_size
    return t_ang


def ang2frac(t_ang, pix_size):
    t_frac = np.asarray(t_ang, dtype="float") / pix_size
    return t_frac


def get_dim(model, shiftmodel="new1.cif"):
    import gemmi

    st = gemmi.read_structure(model)
    model = st[0]
    com = model.calculate_center_of_mass()
    print(com)
    xc = []
    yc = []
    zc = []
    for cra in model.all():
        cra.atom.pos.x -= com.x
        cra.atom.pos.y -= com.y
        cra.atom.pos.z -= com.z
        xc.append(cra.atom.pos.x)
        yc.append(cra.atom.pos.y)
        zc.append(cra.atom.pos.z)
    st.spacegroup_hm = "P 1"
    st.make_mmcif_document().write_file(shiftmodel)
    xc_np = np.asarray(xc)
    yc_np = np.asarray(yc)
    zc_np = np.asarray(zc)
    distances = np.sqrt(np.power(xc_np, 2) + np.power(yc_np, 2) + np.power(zc_np, 2))
    dim1 = 2 + (int(np.max(distances)) + 1) * 2
    return dim1

def rm_zyx2xyz(op):
    assert op.ndim == 2
    assert op.shape[0] == op.shape[1] == 3
    tmp = np.zeros(op.shape, 'float')
    rm = np.zeros(op.shape, 'float')
    tmp[:,0] = op[:,2]
    tmp[:,1] = op[:,1]
    tmp[:,2] = op[:,0]
    rm[0, :] = tmp[2, :]
    rm[1, :] = tmp[1, :]
    rm[2, :] = tmp[0, :] 
    return rm
