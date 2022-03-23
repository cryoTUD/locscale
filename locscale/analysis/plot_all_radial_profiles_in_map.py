#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 25 19:27:12 2022

@author: alok
"""
#%% Introduction
'''
This script is used to plot the various radial profiles present locally in the map

Input: 
    1) atomic model path
    2) mask path
    
Output: 
    1) python dictionary with 10000 center positions and corresponding radial profiles
    
'''

import mrcfile
import gemmi
import os
import numpy as np
import random
from tqdm import tqdm
from locscale.include.emmer.pdb.pdb_tools import find_wilson_cutoff
from locscale.include.emmer.pdb.pdb_to_map import pdb2map
from locscale.include.emmer.pdb.pdb_utils import shake_pdb, set_atomic_bfactors
from locscale.include.emmer.pdb.pdb_tools import get_all_atomic_positions, get_atomic_bfactor_window
from locscale.include.emmer.ndimage.map_utils import extract_window, convert_pdb_to_mrc_position, resample_image
from locscale.include.emmer.ndimage.profile_tools import compute_radial_profile, frequency_array, plot_radial_profile
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
#%%% Inputs
folder = "/mnt/c/Users/abharadwaj1/Downloads/ForUbuntu/faraday_discussions/Tests/test_bfactor_assignment"
rmsd_magnitudes = [0,200,2000]
model_map_paths = {}
model_maps = {}
for rmsd in rmsd_magnitudes:
    model_map_paths[rmsd] = os.path.join(folder, "pdb6y5a_perturbed_{}_pm.mrc".format(rmsd))
    model_maps[rmsd] = mrcfile.open(model_map_paths[rmsd]).data

mask_path = os.path.join(folder, "pdb6y5a_refmac_refined_emd_additional_model_mask.mrc")  ### Ensure that mask doesnt' have unmodelled regions
mask = mrcfile.open(mask_path).data


## Simulation parameters

apix =mrcfile.open(mask_path).voxel_size.tolist()[0]   # angstoerm per pixel
map_shape = mask.shape   ## (voxels, voxels, voxels)
boxsize_real_length = 25 ## angstoerm

boxsize = int(round(boxsize_real_length / apix))
nyquist_freq = apix*2

sample_size = 5000

all_inside_mask = np.asarray(np.where(mask>=0.99)).T.tolist()

random_sample = random.sample(all_inside_mask, sample_size)

bfactor_profiles = {}


for rmsd in model_maps.keys():
    model_map = model_maps[rmsd]
    temp_dictionary = {}
    for center in random_sample:

        center = tuple(center)
        map_window = extract_window(model_map, center, boxsize)
        rp_window = compute_radial_profile(map_window)
        freq = frequency_array(rp_window, apix)
        assert len(freq) == len(rp_window)
    
        temp_dictionary[center] = rp_window
    
    bfactor_profiles[rmsd] = temp_dictionary


#%%

def plot_list_radial_profile(freq,list_of_list_of_profiles,legends=None, font=22,showlegend=True, showPoints=True, alpha=0.1, variation=None, logScale=True):
    import matplotlib.pyplot as plt
    from matplotlib.pyplot import cm
    
    plt.rc('font',size=font)
    
    color_shade = ["red","blue","green"]
    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    ax1.grid(False)
    ax2 = ax1.twiny()
    
    for i,list_of_profiles in enumerate(list_of_list_of_profiles):

        profile_list = np.array(list_of_profiles)
        average_profile = np.einsum("ij->j", profile_list) / len(profile_list)
            
        variation = []
        for col_index in range(profile_list.shape[1]):
            col_extract = profile_list[:,col_index]
            variation.append(col_extract.std())
    
        variation = np.array(variation)
            
        y_max = average_profile + variation
        y_min = average_profile - variation
    
        
            

            
        if logScale:
            ax1.plot(freq**2, np.log(y_max), color_shade[i][0],alpha=1)
            ax1.plot(freq**2, np.log(y_min), color_shade[i][0],alpha=1)
            ax1.fill_between(freq**2,np.log(y_max), np.log(y_min), color=color_shade[i], alpha=alpha)
            
            ax1.legend(["N={}".format(len(profile_list))])
            
                
            ax2.set_xticks(ax1.get_xticks())
            ax2.set_xbound(ax1.get_xbound())
            ax2.set_xticklabels([round(1/np.sqrt(x),1) for x in ax1.get_xticks()])
                
        
            ax1.set_xlabel(r'$1/d^2 [\AA^{-2}]$',fontsize=font)
            ax1.set_ylabel('$ln\mid F \mid $',fontsize=font)
            ax2.set_xlabel('$d [\AA]$',fontsize=font)
        else:
            ax1.plot(freq, y_max, color_shade[i][0],alpha=1)
            ax1.plot(freq, y_min, color_shade[i][0],alpha=1)
            ax1.fill_between(freq,y_max, y_min,color=color_shade[i], alpha=alpha)
            
            ax1.legend(["N={}".format(len(profile_list))])
            
                    
            ax2.set_xticks(ax1.get_xticks())
            ax2.set_xbound(ax1.get_xbound())
            ax2.set_xticklabels([round(1/x,1) for x in ax1.get_xticks()])
            
        
            ax1.set_xlabel('$1/d [\AA^{-1}]$',fontsize=font)
            ax1.set_ylabel('normalised $\mid F \mid $',fontsize=font)
            ax2.set_xlabel('$d [\AA]$',fontsize=font)
        
  
    
    plt.tight_layout()
    return fig

#%%
plot_list_radial_profile(freq, list_of_list_of_profiles=[list(bfactor_profiles[x].values()) for x in [2000]], logScale=True, alpha=0.1)


#%%

# Extract amplitudes at each frequency bins for a given list of profiles

def extract_amplitudes_per_frequency_bins(freq_array, list_of_profiles):
    amplitudes_frequency_bin = {}
    array_of_profiles = np.array(list_of_profiles)
    amplitudes_freq = {}
    for col_index in range(len(freq_array)):
        freq = freq_array[col_index]
        col_extract = array_of_profiles[:,col_index]
        amplitudes_freq[freq] = col_extract
    
    return amplitudes_freq

amplitudes_frequency_bins_rmsd = {}

for rmsd in rmsd_magnitudes:
    amplitudes_frequency_bins_rmsd[rmsd] = extract_amplitudes_per_frequency_bins(freq, list(bfactor_profiles[rmsd].values()))
#%%

def compare_distributions(array1, array2):
    from scipy.stats import pearsonr, ks_2samp
    from scipy.stats import gaussian_kde
    import random
    
    
    array1_local = array1.copy()
    array2_local = array2.copy()
    random.shuffle(array1_local)
    random.shuffle(array2_local)
    kde1 = gaussian_kde(array1_local)
    kde2 = gaussian_kde(array2_local)
    
    bins = np.linspace(min(array1_local.min(),array2_local.min()), max(array1_local.max(),array2_local.max()), 200)
    
    pdf1 = kde1(bins)
    pdf2 = kde2(bins)

    
    r2 = pearsonr(pdf1, pdf2)[0]
    return r2

def get_amplitude_histogram_correlation_curve(rmsd1, rmsd2):
    corr = {}
    for freq_index in range(len(freq)):
        freq_value = freq[freq_index]
        c = compare_distributions(amplitudes_frequency_bins_rmsd[rmsd1][freq[freq_index]], amplitudes_frequency_bins_rmsd[rmsd2][freq[freq_index]])
        corr[freq_value] = c
    
    return corr

amplitude_histogram_correlation = get_amplitude_histogram_correlation_curve(0,2000)
plt.plot(amplitude_histogram_correlation.keys(), amplitude_histogram_correlation.values(),'k--.')


#%%
freq_index = -2
sns.kdeplot(amplitudes_frequency_bins_rmsd[0][freq[freq_index]]),sns.kdeplot(amplitudes_frequency_bins_rmsd[200][freq[freq_index]]),sns.kdeplot(amplitudes_frequency_bins_rmsd[2000][freq[freq_index]]), plt.legend(["0A","2A","20A"]), plt.title("Freq = {}".format(round(1/freq[freq_index],1)))
        
    
    
