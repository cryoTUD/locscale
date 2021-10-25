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
    
    emmap = mrcfile.open(emmap_path).data
    apix = mrcfile.open(emmap_path).voxel_size.tolist()[0]
    rp_emmap = compute_radial_profile(emmap)
    freq = frequency_array(rp_emmap, apix=apix)
    
    bfactor, amp, (fit, z, slope) = estimate_bfactor_through_pwlf(freq, rp_emmap, wilson_cutoff, fsc_cutoff, num_segments=num_segments)
    debye_slope = slope[1]
    print("Debye slope is: ",debye_slope)
    print("Breakpoints and slopes: ",1/np.sqrt(z), slope)
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
    simmap = pdb2map(st, apix=apix, size=size, verbose=False, set_refmac_blur=True)
    
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
        
        plt.plot(freq, fsc_curve,'b')
        plt.plot(freq, np.ones(len(fsc_curve))*metric,'r--')
        
    print("Map quality measured by {} is {}".format(test, round(metric,2)))
    return metric


def obtain_map_quality_multiple(list_of_emmap_path, mask_path, pdb_path):
    map_quality_rscc = {}
    for emmap_path in list_of_emmap_path:
        map_name = emmap_path.split("/")[-1]
        map_quality_rscc[map_name] = map_quality_pdb(emmap_path, mask_path, pdb_path)
    
    print(map_quality_rscc)
    

def plot_fsc_metric_multiple(list_of_emmap_path, mask_path, pdb_path):
    from locscale.include.emmer.ndimage.fsc_util import calculate_fsc_maps
    from locscale.include.emmer.pdb.pdb_utils import set_atomic_bfactors
    from locscale.include.emmer.pdb.pdb_to_map import pdb2map
    from locscale.include.emmer.ndimage.profile_tools import frequency_array
    fsc_quality_multiple = {}
    mask = mrcfile.open(mask_path).data
    apix = mrcfile.open(mask_path).voxel_size.tolist()[0]
    size=mask.shape
    st = gemmi.read_structure(pdb_path)
    st_0 = set_atomic_bfactors(input_gemmi_st=st, b_iso=0)
    simmap = pdb2map(st, apix=apix, size=size, verbose=False, set_refmac_blur=True)
    
    masked_simmap = mask * simmap
    
    fsc_curves = {}
    map_names_legend = []
    fsc_metric = {}
    for emmap_path in list_of_emmap_path:
        from locscale.include.emmer.ndimage.map_utils import resample_image
        map_name = emmap_path.split("/")[-1]
        print(map_name)
        emmap = mrcfile.open(emmap_path).data
        masked_emmap = mask * emmap
        fsc_curves[map_name] = calculate_fsc_maps(masked_emmap, masked_simmap)
        fsc_metric[map_name] = fsc_curves[map_name].mean()
        map_names_legend.append(map_name)
    
    from locscale.include.emmer.ndimage.fsc_util import plot_fscs
    
    freq = frequency_array(fsc_curves[map_name], apix=apix)
    plot_fscs(freq, list(fsc_curves.values()), legends=list(fsc_curves.keys()))
    
    print(fsc_metric)
    
    
    
        
        
        
        