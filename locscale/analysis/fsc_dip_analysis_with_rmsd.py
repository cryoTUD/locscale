#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 25 19:27:12 2022

@author: alok
"""
#%% Introduction
'''
This script is used to find an "Average local radial profile" of a refined atomic model

Input: 
    1) atomic model path
    2) mask path
    3) scattering magnitude = 10 (default)
Output: 
    1) python dictionary with 10000 center positions and corresponding radial profiles
    2) local bfactor correlation of scattered atomic model map and atomic bfactors of scattered model
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
rmsd = [0,200,2000]
model_map_paths = [os.path.join(folder, "pdb6y5a_rmsd_{}_pm_perturbed_using_mask_shifted_4locscale.mrc".format(x)) for x in rmsd]

mask_path = os.path.join(folder, "pdb6y5a_refmac_refined_emd_additional_model_mask.mrc")  ### Ensure that mask doesnt' have unmodelled regions
mask = mrcfile.open(mask_path).data

perturb_magnitude = 10   ## in angstoerms

## Simulation parameters

apix =mrcfile.open(mask_path).voxel_size.tolist()[0]   # angstoerm per pixel
map_shape = mask.shape   ## (voxels, voxels, voxels)
boxsize_real_length = 25 ## angstoerm

boxsize = int(round(boxsize_real_length / apix))
nyquist_freq = apix*2

sample_size = 1000

all_inside_mask = np.asarray(np.where(mask>=0.99)).T.tolist()

random_sample = random.sample(all_inside_mask, sample_size)

bfactor_profiles = {}

##  TO BE CONTINUED!! ADD MODEL MAPS HERE
for center in tqdm(random_sample, desc="Extracting radial profiles"):
    try: 
        center = tuple(center)
        map_window = extract_window(simulated_map, center, boxsize)
        rp_window = compute_radial_profile(map_window)
        freq = frequency_array(rp_window, apix)
        assert len(freq) == len(rp_window)
    
        bfactor_profiles[center] = rp_window
    except:
        continue
    
#%% Plot

profile_list = np.array(list(bfactor_profiles.values()))

plot_radial_profile(freq, list_of_profiles=profile_list)

def extract_local_radial_profiles(input_map, input_mask, apix, boxsize, sample_size=10000, mask_threshold=0.99):
    '''
    Function to find average local radial profile by analysing local profiles at large number of points

    Parameters
    ----------
    input_map : emmap path or emmap
        DESCRIPTION
    input_mask : TYPE
        DESCRIPTION.
    sample_size : TYPE, optional
        DESCRIPTION. The default is 10000.
    boxsize : int, 
        Size of the local window in pixels to analyse the radial profiles

    Returns
    -------
    None.

    '''
    
    from locscale.include.emmer.ndimage.map_utils import extract_window, parse_input
    from locscale.include.emmer.ndimage.profile_tools import compute_radial_profile, frequency_array
    import random
   
    
    emmap = parse_input(input_map)
    mask = parse_input(input_mask)
    all_inside_mask = np.asarray(np.where(mask>=mask_threshold)).T.tolist()
  
    
    bfactor_profiles = {}
    
    num_profiles = 0
    list_of_centers = []
    expected_length_radial_profile = int(boxsize/2)
    print("Extracting profiles...")
    while num_profiles < sample_size:
        try:
            center = random.choice(all_inside_mask)
            if center not in list_of_centers:
                center = tuple(center)
                map_window = extract_window(emmap, center, boxsize)
                rp_window = compute_radial_profile(map_window)
                assert len(rp_window) == expected_length_radial_profile
                list_of_centers.append(center)
                bfactor_profiles[center] = rp_window
                num_profiles += 1
            else:
                continue
        except:
            continue
        
    assert len(list_of_centers) == sample_size
    
    local_profiles = np.array(list(bfactor_profiles.values()))
    freq_array = frequency_array(expected_length_radial_profile, apix)
    
    return (freq_array, local_profiles)

    

    
    


