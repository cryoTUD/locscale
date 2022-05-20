#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 19 16:16:40 2021

@author: alok
"""
import mrcfile
import numpy as np
import pandas as pd
import gemmi

def measure_debye_pwlf(emmap_path, wilson_cutoff, fsc_cutoff, num_segments=3, plot_profile=False, plot_legends=None):
    import mrcfile
    from locscale.include.emmer.ndimage.profile_tools import compute_radial_profile, frequency_array, estimate_bfactor_through_pwlf, plot_radial_profile
    
    if fsc_cutoff > 5:
        print("Resolution too poor to estimate debye slope. Returning zero")
        return 0
    
    emmap = mrcfile.open(emmap_path).data
    apix = mrcfile.open(emmap_path).voxel_size.tolist()[0]
    rp_emmap = compute_radial_profile(emmap)
    freq = frequency_array(rp_emmap, apix=apix)
    
    bfactor, amp, (fit, z, slope) = estimate_bfactor_through_pwlf(freq, rp_emmap, wilson_cutoff, fsc_cutoff, num_segments=num_segments)

    debye_slope = abs(slope[1]-slope[2])

    print("Debye slope is: ",debye_slope)
    print("Breakpoints and slopes: ",1/np.sqrt(z), slope)
    print("Fit quality", round(fit.r_squared(), 2))
    if plot_profile:
        import matplotlib.pyplot as plt
        rp_predict = np.exp(fit.predict(np.copy(freq**2)))
        fig = plot_radial_profile(freq,[rp_emmap, rp_predict])
        if plot_legends is not None:
            plt.legend(plot_legends)
            
    return debye_slope
    

def map_quality_kurtosis(emmap_path, mask_path=None):
    from scipy.stats import kurtosis
    emmap = mrcfile.open(emmap_path).data
    if mask_path is not None:
        mask = mrcfile.open(mask_path).data
        emmap = emmap * mask
    k = kurtosis(emmap.flatten())
    print("Map kurtosis is: {}".format(round(k,2)))
    return k

def map_quality_pdb_multiple(list_of_emmap_path, mask_path, pdb_path):
    from locscale.include.emmer.ndimage.map_tools import compute_real_space_correlation as rscc
    from locscale.include.emmer.ndimage.fsc_util import calculate_fsc_maps
    from locscale.include.emmer.pdb.pdb_utils import set_atomic_bfactors
    from locscale.include.emmer.pdb.pdb_to_map import pdb2map
    from locscale.include.emmer.ndimage.profile_tools import frequency_array
    
    metrics_per_emmap = {}
    mask = mrcfile.open(mask_path).data
    apix = mrcfile.open(mask_path).voxel_size.tolist()[0]
    size=mask.shape
    st = gemmi.read_structure(pdb_path)
    st_0 = set_atomic_bfactors(input_gemmi_st=st, b_iso=0)
    simmap = pdb2map(st_0, apix=apix, size=size, verbose=False, set_refmac_blur=True)
    masked_simmap = simmap * mask
        
    for emmap_path in list_of_emmap_path:
        metric_dictionary = {}
        emmap = mrcfile.open(emmap_path).data

        
        masked_emmap = emmap * mask
        
        
        metric_dictionary['rscc'] = rscc(masked_emmap, masked_simmap)
 
        fsc_curve = calculate_fsc_maps(masked_emmap, masked_simmap)
        freq = frequency_array(fsc_curve, apix=apix)
        metric_dictionary['fsc'] = fsc_curve.mean()
        metrics_per_emmap[emmap_path] = metric_dictionary
        print(emmap_path.split("/")[-1], metrics_per_emmap[emmap_path])

    return metrics_per_emmap

def map_quality_pdb(emmap_path, mask_path, pdb_path, test='rscc'):
    from locscale.include.emmer.ndimage.map_tools import compute_real_space_correlation as rscc
    from locscale.include.emmer.ndimage.fsc_util import calculate_fsc_maps
    from locscale.include.emmer.pdb.pdb_utils import set_atomic_bfactors
    from locscale.include.emmer.pdb.pdb_to_map import pdb2map
    
    emmap = mrcfile.open(emmap_path).data
    mask = mrcfile.open(mask_path).data
    apix = mrcfile.open(emmap_path).voxel_size.tolist()[0]
    size=emmap.shape
    st = gemmi.read_structure(pdb_path)
    st_0 = set_atomic_bfactors(input_gemmi_st=st, b_iso=0)
    simmap = pdb2map(st_0, apix=apix, size=size, verbose=False, set_refmac_blur=True)
    
    masked_emmap = emmap * mask
    masked_simmap = simmap * mask
    
    if test=='rscc':
        metric = rscc(masked_emmap, masked_simmap)
        
    
    if test=='fsc':
        from locscale.include.emmer.ndimage.profile_tools import frequency_array
        import matplotlib.pyplot as plt
        
        fsc_curve = calculate_fsc_maps(masked_emmap, masked_simmap)
        freq = frequency_array(fsc_curve, apix=apix)
        metric = fsc_curve.mean()
        
        
    print("Map quality measured by {} is {}".format(test, round(metric,2)))
    return metric

def local_histogram_analysis(emmap_path, mask_path, fsc_resolution, num_locations=10000, window_size=40):
    import mrcfile
    import random
    from locscale.include.emmer.ndimage.profile_tools import estimate_bfactor_standard, compute_radial_profile, frequency_array, plot_radial_profile
    from locscale.include.emmer.pdb.pdb_tools import find_wilson_cutoff
    from scipy.stats import kurtosis, skew
    from tqdm import tqdm
    
    def get_box(big_volume,center,size):
        return big_volume[center[2]-size//2:center[2]+size//2,center[1]-size//2:center[1]+size//2,center[0]-size//2:center[0]+size//2]

    def distance_from_center_of_box(center_of_window,shape):
        zw,yw,xw = center_of_window
        zc,yc,xc = shape[0]//2, shape[1]//2, shape[2]//2
        
        r = np.sqrt((zc-zw)**2 + (yc-yw)**2 + (xc-xw)**2)
        
        return r
    
    def linear(x,a,b):
        return a * x + b

    def general_quadratic(x,a,b,c):
        return a * x**2 + b*x + c
        
    def r2(y_fit, y_data):
        y_mean = y_data.mean()
        residual_squares = (y_data-y_fit)**2
        variance = (y_data-y_mean)**2
        
        residual_sum_of_squares = residual_squares.sum()
        sum_of_variance = variance.sum()
        
        r_squared = 1 - residual_sum_of_squares/sum_of_variance
        
        return r_squared
    
    def regression(data_input, x_col, y_col,kind="linear"):
        from scipy.optimize import curve_fit
        
        data_unsort = data_input.copy()
        data=data_unsort.sort_values(by=x_col)
        x_data = data[x_col]
        y_data = data[y_col]
        
        if kind == "linear":
            p_opt, p_cov = curve_fit(linear, x_data, y_data)
            a,b = p_opt
            y_fit = linear(x_data, *p_opt)
            r_squared = r2(y_fit, y_data)
            
            return r_squared
        elif kind == "quadratic":
                
            p_opt, p_cov = curve_fit(general_quadratic, x_data, y_data)
            a,b,c = p_opt
            y_fit = general_quadratic(x_data, *p_opt)
            r_squared = r2(y_fit, y_data)
            
            return r_squared
        
        else:
            return None
    

    mask = mrcfile.open(mask_path).data    
    emmap_unmasked = mrcfile.open(emmap_path).data
    apix = mrcfile.open(emmap_path).voxel_size.tolist()[0]
    emmap = emmap_unmasked * mask
    
    wilson_cutoff = find_wilson_cutoff(mask_path=mask_path)
    fsc_cutoff = fsc_resolution
    
    z,y,x = np.where(mask == 1)
    all_points = list(zip(x,y,z))
    random_centers = random.sample(all_points,num_locations)
    
    local_analysis = {}
    
    for center in tqdm(random_centers, desc="Local analysis"):
        try:
            distance_to_center = distance_from_center_of_box(center, emmap.shape)

            window_emmap = get_box(emmap, center, window_size)
            
            ## calculate rp

            rp_emmap = compute_radial_profile(window_emmap)
            freq = frequency_array(rp_emmap, apix)
            
            bfactor_emmap = estimate_bfactor_standard(freq, rp_emmap, wilson_cutoff=wilson_cutoff, fsc_cutoff=fsc_cutoff)
            
            ## histogram metrics
            
 
            skew_emmap = skew(window_emmap.flatten())
            kurtosis_emmap = kurtosis(window_emmap.flatten())
        

            mean_emmap = window_emmap.mean()
            variance_emmap = window_emmap.var()
            
            
            local_analysis[center] = [bfactor_emmap, mean_emmap, variance_emmap, kurtosis_emmap, skew_emmap, tuple(center), distance_to_center]
        except:
            continue
    
    df = pd.DataFrame(data=local_analysis.values(), columns=['bfactor_emmap','mean_emmap','variance_emmap', 'kurtosis_emmap', 'skew_emmap', 'center', 'radius'])
    
    r_squared_skew_kurtosis = regression(df, 'skew_emmap', 'kurtosis_emmap', kind="quadratic")
    r_squared_mean_variance = regression(df, 'mean_emmap', 'variance_emmap', kind="linear")
    
    return df, r_squared_skew_kurtosis, r_squared_mean_variance
          
        
        
        