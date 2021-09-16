# test rotation axis
import sys
import numpy as np
import fcodes_fast as fcodes
import emda.emda_methods as em
from scipy import linalg
from emda.core.quaternions import isRotationMatrix

debug_mode = 0


def map_cov(arr):
    cov = np.zeros(9, dtype="float").reshape(3, 3)
    grid = np.indices(arr.shape)
    for i in range(3):
        for j in range(3):
            cov[-i - 1, -j - 1] = np.sum(
                (grid[i] - np.average(grid[i])) * (grid[j] - np.average(grid[j])) * arr
            ) / np.sum(arr)
    return cov


def apply_rotation(arr, rotmat):
    nx, ny, nz = arr.shape
    rotated_map = fcodes.trilinear_map(rotmat, arr, debug_mode, nx, ny, nz)
    return rotated_map


if __name__ == "__main__":

    #refmap = "/Users/ranganaw/MRC/REFMAC/tutorial_data/reference_map.mrc"
    #targetmap = "/Users/ranganaw/MRC/REFMAC/tutorial_data/emd_10574_fullmap_resampled.mrc"

    refmap, targetmap = sys.argv[1:]

    uc, arr1, origin = em.read_map(refmap)
    uc, arr2, origin = em.read_map(targetmap)
    arr_list = [arr1, arr2]

    nx, ny, nz = arr1.shape
    box_center = (nx//2, ny//2, nz//2)
    eigenvecs = []
    eigenvals = []
    arr_mvd_list = []
    for arr in arr_list:
        com1 = em.center_of_mass_density(arr)
        arr_mvd = em.shift_density(arr, np.subtract(box_center, com1))
        arr_mvd_list.append(arr_mvd)
        map_cov1 = map_cov(arr_mvd * (arr_mvd > np.average(arr_mvd) * 0.4))
        eigvals1, eigvecs1 = linalg.eig(map_cov1)
        eigenvecs.append(eigvecs1)
        eigenvals.append(eigvals1)
        print(eigvals1)
        print(eigvecs1)
    # temporary modification for EMD-10574 map on the reference map
    v1 = eigenvecs[0]
    v2 = eigenvecs[1]
    tmp = np.zeros(v2.shape, 'float')
    tmp[:,0] = v2[:,1]
    tmp[:,1] = -v2[:,0]
    tmp[:,2] = v2[:,2]
    r_ab = tmp @ v1
    #print(isRotationMatrix(r_ab))
    # output PCA applied map
    em.write_mrc(apply_rotation(arr_mvd_list[1], r_ab), 'fitted_map_pca.map', uc)
    