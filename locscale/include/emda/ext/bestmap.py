"""
Author: "Rangana Warshamanage, Garib N. Murshudov"
MRC Laboratory of Molecular Biology
    
This software is released under the
Mozilla Public License, version 2.0; see LICENSE.
"""

# calcuating bestmap - normalized and FSC weighted map

def bestmap(f1, f2, mode=1, kernel_size=5, bin_idx=None, nbin=None, res_arr=None, B=None):
    if mode == 1:
        print("bestmap is calculated using FSC in resolution bins")
        bestmap = bestmap_1dfsc(f1, f2, bin_idx, nbin, res_arr, B)
    elif mode == 2:
        print("bestmap is calculated using local FSC")
        bestmap = bestmap_3dfsc(f1, f2, kernel_size)
    return bestmap


def bestmap_1dfsc(f1, f2, bin_idx, nbin, res_arr=None, B=None):
    import fcodes_fast
    from emda.core import fsc, iotools
    import numpy as np
    import pandas

    cx, cy, cz = f1.shape
    fsc, var_n, var_s, vat_t, _, eo = fsc.halfmaps_fsc_variance(f1, f2, bin_idx, nbin)
    if res_arr is not None:
        data = np.column_stack((res_arr, var_n, var_s, fsc))
    else:
        idx_arr = np.arange(nbin, dtype='int')
        data = np.column_stack((idx_arr, var_n, var_s, fsc))
    columns = ['resol/indx', 'var_n', 'var_s', 'FSC_half']
    df = pandas.DataFrame(data=data, columns=columns)
    iotools.output_to_table(df)
    fsc = 2 * fsc / (1 + fsc)

    if B is None:
        fsc_grid = fcodes_fast.read_into_grid(bin_idx, fsc, nbin, cx, cy, cz)
        fsc_grid_filtered = np.where(fsc_grid < 0.0, 0.0, fsc_grid)
        return np.sqrt(fsc_grid_filtered) * eo
    else:
        k2 = np.exp(-B/res_arr**2/2)
        fac = np.sqrt(k2)*np.sqrt(np.where(fsc<0, 0, fsc))/(1+(k2-1)*fsc)
        fac_grid = fcodes_fast.read_into_grid(bin_idx, fac, nbin, cx, cy, cz)
        return fac_grid * eo 


def bestmap_3dfsc(f1, f2, kernel_size=5):
    import numpy as np
    from emda.core.restools import create_soft_edged_kernel_pxl
    from emda.ext.fouriersp_local import get_3dnoise_variance, get_3dtotal_variance

    kernel = create_soft_edged_kernel_pxl(kernel_size)
    # calculate 3D fourier correlation
    noise = get_3dnoise_variance(f1, f2, kernel)
    total = get_3dtotal_variance(f1, f2, kernel)
    total_pos = np.where(total <= 0.0, 0.1, total)
    fsc = (total - noise) / total_pos
    fsc_grid_filtered = np.where((total - noise) <= 0.0, 0.0, fsc)
    eo = (f1 + f2) / (2 * np.sqrt(total_pos))
    eo = np.where(total <= 0.0, 0.0, eo)
    return np.sqrt(fsc_grid_filtered) * eo
