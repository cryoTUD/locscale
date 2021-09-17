# symmetrize map by operators
import numpy as np
import emda.emda_methods as em
import fcodes_fast
from emda.core.restools import get_resolution_array
from emda.ext.mapfit.utils import double_the_axes
from numpy.fft import fftn, ifftn, fftshift, ifftshift
from emda.ext.sym.GenerateOperators_v9_ky4 import from_two_axes_to_group_v2, operators_from_symbol


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
    return frs

def apply_op2(f1, op, bin_idx, nbin):
    nz, ny, nx = f1.shape 
    frs = fcodes_fast.trilinear2(f1,bin_idx,op,nbin,0,1,nz,ny,nx)[:,:,:,0]
    return frs

def get_matrices(foldlist, axlist):
    if len(foldlist) == 2:
        order1, order2 = foldlist[0], foldlist[1]
        ax1, ax2 = axlist[0], axlist[1]
        ax1, ax2 = ax1[::-1], ax2[::-1]
        pg, _, _, ops = from_two_axes_to_group_v2(axis1=ax1, order1=order1, axis2=ax2, order2=order2)
        print('PG= ', pg)
    if len(foldlist) == 1:
        order1 = foldlist[0]
        ax1 = axlist[0]
        ax1 = ax1[::-1] # reverse the axis
        pg, _, _, ops = from_two_axes_to_group_v2(axis1=ax1, order1=order1)
        print('PG= ', pg)
    return ops

def rebox_map(arr1):
    nx, ny, nz = arr1.shape
    dx = int(nx / 4)
    dy = int(ny / 4)
    dz = int(nz / 4)
    reboxed_map = arr1[dx : dx + nx//2, dy : dy + ny//2, dz : dz + nz//2]
    return reboxed_map


def symmetrize_map_known_pg(imap, pg, outmapname=None):
    if outmapname is None:
        outmapname = "emda_sym_averaged_map.mrc"
    _, _, ops = operators_from_symbol(pg)
    uc, arr, orig = em.get_data(imap)
    arr2 = double_the_axes(arr)
    f1 = fftshift(fftn(fftshift(arr2)))
    nbin, res_arr, bin_idx = get_resolution_array(uc, f1)
    frs_sum = f1
    print("Symmetrising map...")
    for op in ops[1:]:
        frs = apply_op(f1, op, bin_idx, nbin)
        frs_sum += frs
    avg_f = frs_sum / len(ops)
    avgmap = ifftshift(np.real(ifftn(ifftshift(avg_f))))
    avgmap = rebox_map(avgmap)
    em.write_mrc(avgmap, outmapname, uc, orig)
    return [avgmap, uc, orig]


def symmetrize_map_using_ops(imap, ops, outmapname=None):
    if outmapname is None:
        outmapname = "emda_sym_averaged_map.mrc"
    uc, arr, orig = em.get_data(imap)
    arr2 = double_the_axes(arr)
    f1 = fftshift(fftn(fftshift(arr2)))
    nbin, res_arr, bin_idx = get_resolution_array(uc, f1)
    frs_sum = f1
    print("Symmetrising map...")
    for op in ops:
        frs = apply_op(f1, op, bin_idx, nbin)
        frs_sum += frs
    avg_f = frs_sum / (len(ops) + 1)
    avgmap = ifftshift(np.real(ifftn(ifftshift(avg_f))))
    avgmap = rebox_map(avgmap)
    em.write_mrc(avgmap, outmapname, uc, orig)
    return [avgmap, uc, orig]




