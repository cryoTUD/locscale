# symmetrize map by operators
import numpy as np
import fcodes_fast
from numpy.fft import fftn, ifftn, fftshift, ifftshift
from locscale.include.symmetry_emda.GenerateOperators_v9_ky4 import operators_from_symbol
from locscale.include.symmetry_emda.apply_rotation_matrix import *
"""
Original authors

Author: "Rangana Warshamanage, Garib N. Murshudov"
MRC Laboratory of Molecular Biology

https://gitlab.com/ccpem/emda/-/tree/master/
EMDA version 1.1.3.post6
"""

def get_resolution_array(uc, hf1):
    debug_mode = 0
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

def double_the_axes(arr1):
    nx, ny, nz = arr1.shape
    big_arr1 = np.zeros((2 * nx, 2 * ny, 2 * nz), dtype="float")
    dx = int(nx / 2)
    dy = int(ny / 2)
    dz = int(nz / 2)
    big_arr1[dx : dx + nx, dy : dy + ny, dz : dz + nz] = arr1
    return big_arr1


def apply_op(f1, op, bin_idx, nbin):
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
    nz, ny, nx = f1.shape 
    frs = fcodes_fast.trilinear2(f1,bin_idx,rm,nbin,0,1,nz,ny,nx)[:,:,:,0]

    #frs = trilinear_interpolation(f1, rm)  # CODE TBC
    #frs = trilinear_interpolation_numpy(f1, rm) # TBC
    return frs

def rebox_map(arr1):
    nx, ny, nz = arr1.shape
    dx = int(nx / 4)
    dy = int(ny / 4)
    dz = int(nz / 4)
    reboxed_map = arr1[dx : dx + nx//2, dy : dy + ny//2, dz : dz + nz//2]
    return reboxed_map


def symmetrize_map_known_pg(emmap, apix, pg):
    print("===== Symmetrize Map =====")
    print("Credits: Rangana Warshamanage, Garib N. Murshudov")
    print("EMDA version 1.1.3.post6")
    print("https://gitlab.com/ccpem/emda/-/tree/master/")
    print("==========================")
    

    _, _, ops = operators_from_symbol(pg)
    #uc, arr, orig = em.get_data(imap)
    unitcell = np.array([emmap.shape[0]*apix, emmap.shape[1]*apix, emmap.shape[2]*apix, 90, 90, 90])
    arr2 = double_the_axes(emmap)
    print("Double the axes: {}".format(arr2.shape))
    f1 = fftshift(fftn(fftshift(arr2)))
    nbin, res_arr, bin_idx = get_resolution_array(unitcell, f1)
    frs_sum = f1
    
    for op in ops:
        frs = apply_op(f1, op,bin_idx,nbin)
        frs_sum += frs
    avg_f = frs_sum / len(ops)
    avgmap = ifftshift(np.real(ifftn(ifftshift(avg_f))))
    avgmap = rebox_map(avgmap)
    
    return avgmap


def symmetrize_map_emda(emmap_path, pg):
    from locscale.include.emmer.ndimage.map_utils import load_map
    emmap,apix = load_map(emmap_path)
    symmetry_average_map = symmetrize_map_known_pg(emmap, apix, pg)
    return symmetry_average_map
