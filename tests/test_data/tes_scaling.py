#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Sep 12 19:02:01 2021

@author: alok
"""

from myheaders import *
from emmer.ndimage.profile_tools import estimate_bfactor_through_pwlf, frequency_array
from pseudomodel_headers import number_of_segments
from compute import compute_radial_profile, compute_scale_factors, set_radial_profile

emmap = mrcfile.open("emd5778_unfiltered.mrc").data
mask = mrcfile.open("emd5778_mask.mrc").data
modmap = mrcfile.open("model_reference.mrc").data
locscale = mrcfile.open("loc_scale.mrc").data

all_inside_mask = np.asarray(np.where(mask>=1)).T.tolist()

random_center = random.sample(all_inside_mask,1)[0]
size=40

emmap_window = extract_window(emmap, random_center, size=size)
modmap_window = extract_window(modmap, random_center, size=size)


rp_emmap = compute_radial_profile(emmap_window)
rp_modmap, radii = compute_radial_profile(modmap_window,return_indices=True)

scale_factor_arguments = {}
scale_factor_arguments['wilson'] = 8.5
scale_factor_arguments['high_freq'] = 4.91
scale_factor_arguments['fsc_cutoff'] = 2.5
scale_factor_arguments['smooth'] = 0.3

scale_factors_old = compute_scale_factors(rp_emmap, rp_modmap, apix=1.2156, scale_factor_arguments=scale_factor_arguments, use_theoretical_profile=False)

scale_factors_new, report = compute_scale_factors(rp_emmap, rp_modmap, apix=1.2156, scale_factor_arguments=scale_factor_arguments, use_theoretical_profile=True, check_scaling=True)

scaled_window_old = set_radial_profile(emmap_window, scale_factors_old, radii=radii)

scaled_window_new = set_radial_profile(emmap_window, scale_factors_new, radii=radii)

rp_scaled_old = compute_radial_profile(scaled_window_old)
rp_scaled_new = compute_radial_profile(scaled_window_new)

freq = frequency_array(rp_emmap, apix=1.2156)
plot_radial_profile(freq, [rp_emmap,rp_modmap, rp_scaled_old, rp_scaled_new, report['scaled_reference_profile'], report['scaled_theoretical_amplitude']])

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

with open("test_scaling_data_3.pickle","wb") as file:
    pickle.dump(test_scaling_data,file)






