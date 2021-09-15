"""
Author: "Rangana Warshamanage, Garib N. Murshudov"
MRC Laboratory of Molecular Biology
    
This software is released under the
Mozilla Public License, version 2.0; see LICENSE.
"""

from __future__ import absolute_import, division, print_function, unicode_literals
from emda.config import *
import numpy as np
from emda import core, ext
import fcodes_fast

# calculate scale between two maps in resolution bins
def scale_onemap2another(f1, f2, bin_idx, res_arr):

    assert f1.shape == f2.shape == bin_idx.shape
    nx, ny, nz = bin_idx.shape
    scale_np = np.zeros(len(res_arr), dtype='float')
    s_grid = fcodes_fast.read_into_grid(bin_idx,
                                1.0/res_arr,
                                len(res_arr),
                                nx,ny,nz)
    print('resolution       scale')
    for i, res in enumerate(res_arr):
        slope = estimate_scale(f1[bin_idx==i], 
                                f2[bin_idx==i], 
                                s_grid[bin_idx==i])
        #scale_np[i] = params[0]
        scale_np[i] = slope
        print("{:8.4f} {:6.2f}".format(res, slope))   
    return scale_np

def estimate_scale(f1, f2, s):
    from emda.ext.mapfit import curve_fit_3
    from scipy import stats
    x0 = np.array([1., 10.])
    if f1.ndim > 1:
        f1 = f1.flatten()
    if f2.ndim > 1:
        f1 = f2.flatten()
    if s.ndim > 1:
        f1 = s.flatten()            
    s = (s**2)/4
    #params = curve_fit_3.lsq(f1, f2, s, x0)
    # just scaling
    slope, intercept,_,_,_ = stats.linregress(np.real(f1*f1.conj()), 
                                              np.real(f2*f2.conj()))
    return slope

def scale_twomaps_by_power(f1, f2, bin_idx=None, uc=None, res_arr=None):
    from emda.ext.mapfit import mapaverage
    assert f1.shape == f2.shape == bin_idx.shape
    nx, ny, nz = f1.shape
    if bin_idx is None:
        nbin,res_arr,bin_idx = core.restools.get_resolution_array(uc,f1)
    else:
        nbin = np.max(bin_idx) + 1
    # find how far the signal is
    bin_stats = core.fsc.anytwomaps_fsc_covariance(f1=f1,
                                                        f2=f2,
                                                        bin_idx=bin_idx,
                                                        nbin=nbin)
    f1f2_fsc, f1f2_covar = bin_stats[0], bin_stats[1]
    mask = (mapaverage.set_array(f1f2_covar, 0.1) > 0.0).astype(int)
    inverse_mask = (mask < 1).astype(int)
    power_1 = fcodes_fast.calc_power_spectrum(f1,bin_idx,nbin,debug_mode,nx,ny,nz)
    power_2 = fcodes_fast.calc_power_spectrum(f2,bin_idx,nbin,debug_mode,nx,ny,nz)
    scale_np = power_2/power_1
    scale_np = scale_np * mask + inverse_mask.astype(float)
    scale_np = 1.0 / scale_np
    core.plotter.plot_nlines_log(res_arr,
                        [power_1,power_2,power_2*scale_np],
                        ["power1","power2","scaled2to1"],
                        'log_totalvariances.eps')
    for i, res in enumerate(res_arr):
        print("{:8.4f} {:6.2f}".format(res, scale_np[i]))
    return scale_np

def transfer_power(bin_idx,res_arr,scale):
    nx, ny, nz = bin_idx.shape
    scale_grid = fcodes_fast.read_into_grid(bin_idx,
                                scale,
                                len(res_arr),
                                nx,ny,nz)
    return scale_grid
