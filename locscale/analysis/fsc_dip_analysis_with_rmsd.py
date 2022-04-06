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
from locscale.include.emmer.ndimage.map_tools import get_atomic_model_mask
from locscale.include.emmer.ndimage.profile_tools import compute_radial_profile, frequency_array, plot_radial_profile
from locscale.include.emmer.ndimage.fsc_util import calculate_fsc_maps
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
#%%% Inputs
input_folder = "/mnt/c/Users/abharadwaj1/Downloads/ForUbuntu/faraday_discussions/Calculations/fsc_dip_analysis"
output_folder = input_folder

pdb_path = os.path.join(input_folder, "pdb6y5a_additional_refined.pdb")
sample_map = os.path.join(input_folder, "refined","locscale_additional_refined_perturbed_strict_masking_no_overlap_rmsd_0_A.mrc")
apix = mrcfile.open(sample_map).voxel_size.tolist()[0]

rmsd_magnitudes = [0,1,2,5,10,15,20]
atomic_model_mask_path = get_atomic_model_mask(emmap_path=sample_map, pdb_path=pdb_path, 
                                               dilation_radius=3, softening_parameter=5)
atomic_model_mask = mrcfile.open(atomic_model_mask_path).data

locscale_blurred_maps = {}
locscale_normal_maps = {}

locscale_prefix_blur = "locscale_additional_refined_perturbed_strict_masking_blurred200_no_overlap_rmsd_{}_A.mrc"
locscale_prefix_refined = "locscale_additional_refined_perturbed_strict_masking_no_overlap_rmsd_{}_A.mrc"

for rmsd in rmsd_magnitudes:
    locscale_blurred_maps[rmsd] = mrcfile.open(os.path.join(input_folder, "blur200", locscale_prefix_blur.format(rmsd))).data
    locscale_normal_maps[rmsd] = mrcfile.open(os.path.join(input_folder, "refined", locscale_prefix_refined.format(rmsd))).data
    


fsc_arrays_blurred = {}
fsc_arrays_normal = {}

for rmsd in rmsd_magnitudes:
    input_map1 = locscale_blurred_maps[0] * atomic_model_mask
    input_map2 = locscale_blurred_maps[rmsd] * atomic_model_mask
    fsc = calculate_fsc_maps(input_map_1=input_map1, input_map_2=input_map2)
    freq = frequency_array(fsc, apix=apix)
    fsc_arrays_blurred[rmsd] = [freq,fsc]

for rmsd in rmsd_magnitudes:
    input_map1 = locscale_normal_maps[0] * atomic_model_mask
    input_map2 = locscale_normal_maps[rmsd] * atomic_model_mask
    fsc2 = calculate_fsc_maps(input_map_1=input_map1, input_map_2=input_map2)
    freq = frequency_array(fsc2, apix=apix)
    fsc_arrays_normal[rmsd] = [freq,fsc]
    
 

   
#%%    
from locscale.include.emmer.ndimage.fsc_util import plot_multiple_fsc
from matplotlib.pyplot import cm


map_pair_tuples = [(locscale_blurred_maps[0], locscale_blurred_maps[rmsd]) for rmsd in rmsd_magnitudes]
map_pair_tuples += [(locscale_normal_maps[0], locscale_normal_maps[rmsd]) for rmsd in rmsd_magnitudes]

legends = ["Blurred: RMSD 0 and {}A".format(rmsd) for rmsd in [0]]
legends += ["Normal: RMSD 0 and {}A".format(rmsd) for rmsd in rmsd_magnitudes]
color_blur = cm.Blues(np.linspace(0,1,len(rmsd_magnitudes)))
color_normal = cm.Reds(np.linspace(0,1,len(rmsd_magnitudes)))
colors = np.concatenate((color_blur, color_normal), axis=0)
plot_multiple_fsc(map_pair_tuples, common_mask_path=atomic_model_mask_path, legend=legends, colors=colors)

#%%
from locscale.utils.plot_utils import pretty_lineplot_multiple_fsc_curves

fsc_dip_blurred = os.path.join(output_folder, "fsc_curves_dip_blurred200_nolegend_with_top_x.eps")
legends=["{} $\AA$".format(rmsd) for rmsd in rmsd_magnitudes]
pretty_lineplot_multiple_fsc_curves(fsc_arrays_blurred, filename=fsc_dip_blurred, fontscale=3, linewidth=2, legends=None)


