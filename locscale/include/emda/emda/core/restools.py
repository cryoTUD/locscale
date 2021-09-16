"""
Author: "Rangana Warshamanage, Garib N. Murshudov"
MRC Laboratory of Molecular Biology
    
This software is released under the
Mozilla Public License, version 2.0; see LICENSE.
"""

from __future__ import absolute_import, division, print_function, unicode_literals

import numpy as np
import fcodes_fast
from emda.config import *


def test():
    print("restools test ... Passed")


def get_resolution_array(uc, hf1):

    nx, ny, nz = hf1.shape
    maxbin = np.amax(np.array([nx // 2, ny // 2, nz // 2]))
    if nx == ny == nz:
        nbin, res_arr, bin_idx, s_grid = fcodes_fast.resol_grid_em(
            uc, debug_mode, maxbin, nx, ny, nz
        )
    else:
        nbin, res_arr, bin_idx, sgrid = fcodes_fast.resolution_grid(
            uc, debug_mode, maxbin, nx, ny, nz
        )
    return nbin, res_arr[:nbin], bin_idx


def get_resArr(unit_cell, nx):
    # Generate resolution array
    a, b, c = unit_cell[:3]
    narr = np.arange(nx // 2)
    narr[0] = 1.0
    fResArr = a / narr
    fResArr[0] = 5000
    return fResArr


def create_kernel(fResArr, smax):
    # Create kernel. smax is resolution to which the kernel
    # is defined.
    dist = np.sqrt((fResArr - smax) ** 2)
    cbin = np.argmin(dist)
    print("cbin", cbin)
    box_size = 2 * cbin + 1
    box_radius = cbin + 1
    # Creating a sphere mask (binary mask)
    center = [box_radius - 1, box_radius - 1, box_radius - 1]
    print("boxsize: ", box_size, "boxradius: ", box_radius, "center:", center)
    radius = box_radius
    X, Y, Z = np.ogrid[:box_size, :box_size, :box_size]
    dist_from_center = np.sqrt(
        (X - center[0]) ** 2 + (Y - center[1]) ** 2 + (Z - center[2]) ** 2
    )
    mask = dist_from_center <= radius
    return mask


def create_binary_kernel(r1):
    from math import sqrt

    boxsize = 2 * r1 + 1
    kern_sphere = np.zeros(shape=(boxsize, boxsize, boxsize), dtype="float")
    kx = ky = kz = boxsize
    center = boxsize // 2
    #print("center: ", center)
    r1 = center
    for i in range(kx):
        for j in range(ky):
            for k in range(kz):
                dist = sqrt((i - center) ** 2 + (j - center) ** 2 + (k - center) ** 2)
                if dist < r1:
                    kern_sphere[i, j, k] = 1
    #kern_sphere = kern_sphere / np.sum(kern_sphere)
    return kern_sphere


def create_soft_edged_kernel(fResArr, smax):
    # Create soft-edged-kernel. smax is resolution to which the kernel size
    # is defined
    from math import sqrt, cos

    mask = create_kernel(fResArr, smax)
    kern_sphere = mask / np.sum(mask)  # Normalized mask
    print("Number of data points inside kernel: ", np.count_nonzero(mask))
    # Creating soft-edged mask - Naive approach
    kern_sphere_soft = np.zeros(shape=(kern_sphere.shape), dtype="float")
    r1 = kern_sphere.shape[0] // 2 + 1
    r0 = r1 - int(round(r1 * 0.3))
    print("r1: ", r1, "r0: ", r0)
    kx = kern_sphere.shape[0]
    ky = kern_sphere.shape[1]
    kz = kern_sphere.shape[2]
    center = [r1 - 1, r1 - 1, r1 - 1]
    for i in range(kx):
        for j in range(ky):
            for k in range(kz):
                dist = sqrt(
                    (i - center[0]) ** 2 + (j - center[0]) ** 2 + (k - center[0]) ** 2
                )
                if dist <= r1:
                    if dist < r0:
                        kern_sphere_soft[i, j, k] = 1
                    else:
                        kern_sphere_soft[i, j, k] = (
                            (1 + cos(np.pi * (dist - r0) / (r1 - r0)))
                        ) / 2.0
    kern_sphere_soft = kern_sphere_soft / np.sum(kern_sphere_soft)
    return kern_sphere_soft


def create_soft_edged_kernel_pxl(r1):
    # Create soft-edged-kernel. r1 is the radius of kernel in pixels
    from math import sqrt, cos

    if r1 < 3:
        boxsize = 5
    else:
        boxsize = 2 * r1 + 1
    kern_sphere_soft = np.zeros(shape=(boxsize, boxsize, boxsize), dtype="float")
    kx = ky = kz = boxsize
    center = boxsize // 2
    #print("center: ", center)
    r1 = center
    r0 = r1 - 2
    #print("r1: ", r1, "r0: ", r0)
    for i in range(kx):
        for j in range(ky):
            for k in range(kz):
                dist = sqrt((i - center) ** 2 + (j - center) ** 2 + (k - center) ** 2)
                if dist < r1:
                    if dist < r0:
                        kern_sphere_soft[i, j, k] = 1
                    else:
                        kern_sphere_soft[i, j, k] = (
                            (1 + cos(np.pi * (dist - r0) / (r1 - r0)))
                        ) / 2.0
    kern_sphere_soft = kern_sphere_soft / np.sum(kern_sphere_soft)
    # contour_nplot2(kern_sphere_soft)
    return kern_sphere_soft


def softedgekernel_5x5():
    from math import sqrt, cos

    boxsize = 5
    kern_sphere_soft = np.zeros(shape=(boxsize, boxsize, boxsize), dtype="float")
    kx = ky = kz = boxsize
    center = boxsize // 2
    r1 = center
    r0 = r1 - 2
    for i in range(kx):
        for j in range(ky):
            for k in range(kz):
                dist = sqrt((i - center) ** 2 + (j - center) ** 2 + (k - center) ** 2)
                if dist < r1:
                    if dist <= r0:
                        kern_sphere_soft[i, j, k] = 1
                    else:
                        kern_sphere_soft[i, j, k] = (
                            (1 + cos(np.pi * (dist - r0) / (r1 - r0)))
                        ) / 2.0
    kern_sphere_soft = kern_sphere_soft / np.sum(kern_sphere_soft)
    return kern_sphere_soft


def remove_edge(fResArr, smax):
    # Remove everything outside smax-resolution.
    dist = np.sqrt((fResArr - smax) ** 2)
    cbin = np.argmin(dist)
    box_radius = cbin + 1
    box_size = cbin * 2 + 1
    # Creating a sphere mask
    center = [box_radius, box_radius, box_radius]
    print("boxsize: ", box_size, "boxradius: ", box_radius, "center:", center)
    radius = box_radius
    X, Y, Z = np.ogrid[:box_size, :box_size, :box_size]
    dist_from_center = np.sqrt(
        (X - center[0]) ** 2 + (Y - center[1]) ** 2 + (Z - center[2]) ** 2
    )
    mask = dist_from_center <= radius
    return mask


def cut_resolution(f, bin_idx, res_arr, cbin):
    # Making data for map fitting
    nx, ny, nz = f.shape
    # fcodes.cutmap cuts f according to resolution defined by smax
    # and output same size map as f but padded with zeros outside.
    fout = fcodes_fast.cutmap(f, bin_idx, cbin, 0, len(res_arr), nx, ny, nz)
    return fout


def cut_resolution_for_linefit(f, bin_idx, res_arr, smax):
    # Making data for map fitting
    nx, ny, nz = f.shape
    cbin = cx = smax
    dx = int((nx - 2 * cx) / 2)
    dy = int((ny - 2 * cx) / 2)
    dz = int((nz - 2 * cx) / 2)
    cBIdx = bin_idx[dx : dx + 2 * cx, dy : dy + 2 * cx, dz : dz + 2 * cx]
    fout = fcodes_fast.cutmap(f, bin_idx, cbin, 0, len(res_arr), nx, ny, nz)[
        dx : dx + 2 * cx, dy : dy + 2 * cx, dz : dz + 2 * cx
    ]
    return fout, cBIdx, cbin


def cut_resolution_nmaps(f_list, bin_idx, res_arr, smax):
    # Making data for map fitting
    f_arr = np.asarray(f_list, dtype="complex")
    nx, ny, nz = f_list[0].shape
    cbin = cx = smax
    print("fit resolution:", res_arr[cbin])
    dx = int((nx - 2 * cx) / 2)
    dy = int((ny - 2 * cx) / 2)
    dz = int((nz - 2 * cx) / 2)
    cBIdx = bin_idx[dx : dx + 2 * cx, dy : dy + 2 * cx, dz : dz + 2 * cx]
    fout = fcodes_fast.cutmap_arr(
        f_arr, bin_idx, cbin, 0, len(res_arr), nx, ny, nz, len(f_list)
    )[:, dx : dx + 2 * cx, dy : dy + 2 * cx, dz : dz + 2 * cx]
    return fout, cBIdx, cbin
