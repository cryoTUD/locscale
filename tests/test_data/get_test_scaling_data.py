#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Sep 12 19:02:01 2021

@author: alok
"""

import mrcfile
import gemmi
from locscale.include.emmer.ndimage.profile_tools import estimate_bfactor_through_pwlf, frequency_array, plot_radial_profile, number_of_segments
from locscale.utils.file_tools import get_locscale_path
from locscale.utils.scaling_tools import compute_scale_factors, set_radial_profile, compute_radial_profile_proper
from locscale.include.emmer.ndimage.map_utils import get_all_voxels_inside_mask, extract_window
from locscale.include.emmer.ndimage.map_tools import compute_radial_profile_simple, set_radial_profile_to_volume
from locscale.include.emmer.pdb.pdb_to_map import pdb2map
from locscale.include.confidenceMapUtil import FDRutil
import os
import random
import numpy as np
import pickle

locscale_path = get_locscale_path()["locscale"]

emmap_path = os.path.join(locscale_path, "tests","test_data","emd5778_map_chainA.mrc")
pdb_path = os.path.join(locscale_path, "tests","test_data","pdb3j5p_refined_chainA.pdb")
mask_path = os.path.join(locscale_path, "tests","test_data","emd5778_mask_chainA.mrc")

MB_locscale_path = os.path.join(locscale_path, "tests","test_data","reference_mb_locscale.mrc")
MF_locscale_path = os.path.join(locscale_path, "tests","test_data","reference_mf_locscale.mrc")

mask = mrcfile.open(mask_path).data

all_inside_mask = get_all_voxels_inside_mask(mask, mask_threshold=0.99)
random_center = random.choice(all_inside_mask)
window_size=40
apix = mrcfile.open(emmap_path).voxel_size.tolist()[0]
emmap = mrcfile.open(emmap_path).data

modmap = pdb2map(input_pdb=pdb_path, apix=apix, size=emmap.shape)
emmap_window = extract_window(emmap, random_center, size=window_size)
modmap_window = extract_window(modmap, random_center, size=window_size)
frequency_map_window = FDRutil.calculate_frequency_map(np.zeros((window_size,window_size,window_size)))


rp_emmap,_ = compute_radial_profile_proper(emmap_window, frequency_map_window)
rp_modmap, frequencies_map = compute_radial_profile_proper(modmap_window, frequency_map_window)

scale_factor_arguments = {}
scale_factor_arguments['wilson'] = 9.69
scale_factor_arguments['nyquist'] = 2.5
scale_factor_arguments['fsc_cutoff'] = 3.4
scale_factor_arguments['smooth'] = 0.3
scale_factor_arguments['boost_secondary_structure'] = 2
scale_factor_arguments['no_reference'] = False

scale_factors_old,_,_ = compute_scale_factors(rp_emmap, rp_modmap, apix=1.2156,  scale_factor_arguments=scale_factor_arguments, use_theoretical_profile=False)
scale_factors_new,_,_,report = compute_scale_factors(rp_emmap, rp_modmap, apix=1.2156, scale_factor_arguments=scale_factor_arguments, use_theoretical_profile=True, check_scaling=True)

scaled_window_old,_ = set_radial_profile(emmap_window, scale_factors_old, frequencies_map, frequency_map_window, emmap_window.shape)
scaled_window_new,_ = set_radial_profile(emmap_window, scale_factors_new, frequencies_map, frequency_map_window, emmap_window.shape)

rp_scaled_old,_ = compute_radial_profile_proper(scaled_window_old, frequency_map_window)
rp_scaled_new,_ = compute_radial_profile_proper(scaled_window_new, frequency_map_window)

freq = frequency_array(rp_emmap, apix=1.2156)
plot_radial_profile(freq, [rp_emmap, rp_modmap, rp_scaled_old, rp_scaled_new],legends=["RP emmap","RP modmap","Scaled classic","Scaled using deviations"], showPoints=False)

test_scaling_data = {}
test_scaling_data['emmap_window'] = emmap_window
test_scaling_data['modmap_window'] = modmap_window
test_scaling_data['scaled_window_old'] = scaled_window_old
test_scaling_data['scaled_window_new'] = scaled_window_new
test_scaling_data['rp_emmap'] = rp_emmap
test_scaling_data['rp_modmap'] = rp_modmap
test_scaling_data['rp_scaled_old'] = rp_scaled_old
test_scaling_data['rp_scaled_new'] = rp_scaled_new
test_scaling_data['report'] = report
test_scaling_data['center'] = random_center

output_pickle_file = os.path.join(locscale_path,"tests","test_data","test_scaling_data_3.pickle")
with open(output_pickle_file,"wb") as file:
    pickle.dump(test_scaling_data,file)






