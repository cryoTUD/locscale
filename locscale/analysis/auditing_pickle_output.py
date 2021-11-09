#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Sep 11 11:17:21 2021

@author: alok
"""
import pickle
import numpy as np
import matplotlib.pyplot as plt
from locscale.include.emmer.ndimage.profile_tools import plot_radial_profile, compute_radial_profile, resample_1d
from locscale.include.emmer.ndimage.map_tools import compute_real_space_correlation as rsc
import os

def get_box(big_volume,center,size):
    return big_volume[center[2]-size//2:center[2]+size//2,center[1]-size//2:center[1]+size//2,center[0]-size//2:center[0]+size//2]

folder = "/mnt/c/Users/abharadwaj1/Downloads/ForUbuntu/LocScale/tests/test_low_frequency_dip"

'''
What's inside the pickle output? For 1% of windows following results are saved:
temporary_dictionary = {}
temporary_dictionary['em_profile'] = em_profile
temporary_dictionary['input_ref_profile'] = ref_profile
temporary_dictionary['freq'] = freq
temporary_dictionary['theoretical_amplitude'] = theoretical_profile_tuple[1]
temporary_dictionary['scaled_theoretical_amplitude'] = scaled_theoretical_amplitude
temporary_dictionary['scaled_reference_profile'] = scaled_reference_profile
temporary_dictionary['scaling_condition'] = [wilson_cutoff, fsc_cutoff]
temporary_dictionary['merging_condition'] = [smooth, wilson_cutoff]
temporary_dictionary['scale_factor'] = scale_factor
'''

pickle_output = os.path.join(folder,"profiles_audit_using_magnification_1p5.pickle")

import mrcfile
model_map = mrcfile.open("/mnt/c/Users/abharadwaj1/Downloads/ForUbuntu/LocScale/tests/test_low_frequency_dip/pdb3j5p_refined_cropped_shifted_refmac_refined_4locscale.mrc").data

with open(pickle_output,"rb") as audit_file:
    audit_scaling = pickle.load(audit_file)

import random
random_positions = list(audit_scaling.keys())    
sample = random.sample(random_positions, 100)
rscc = {
        'interpolated':[],
        'deviation':[]}
for key in sample:
    model_map_window = get_box(model_map, key, size=24)
    
    freq = audit_scaling[key]['freq']
    rp_model_map = compute_radial_profile(model_map_window)
    em_profile = audit_scaling[key]['em_profile']
    ref_profile = audit_scaling[key]['input_ref_profile']
    theoretical_profile = audit_scaling[key]['theoretical_amplitude']
    scaled_theoretical = audit_scaling[key]['scaled_theoretical_amplitude']
    merged_profile = audit_scaling[key]['scaled_reference_profile']
    deviated_profile = audit_scaling[key]['deviated_reference_profile']
    
    _,resampled_model_map = resample_1d(freq, rp_model_map, num=100, xlims=[1/9.6, 1/3.4])
    _,resampled_merged = resample_1d(freq, merged_profile, num=100, xlims=[1/9.6, 1/3.4])
    _,resampled_deviated = resample_1d(freq, deviated_profile, num=100, xlims=[1/9.6, 1/3.4])
    
    rscc['interpolated'].append(rsc(resampled_model_map, resampled_merged))
    rscc['deviation'].append(rsc(resampled_model_map, resampled_deviated))
    '''
    print(audit_scaling[key]['scaling_condition'])
    
    plot_radial_profile(freq,[ref_profile],
                        legends=['ref_profile']);
    
    plot_radial_profile(freq,[scaled_theoretical],
                        legends=['scaled theoretical']);
    
    plot_radial_profile(freq,[ref_profile, scaled_theoretical],
                        legends=['ref_profile','scaled theoretical']);
    
    plot_radial_profile(freq,[ref_profile, merged_profile, scaled_theoretical, deviated_profile],
                        legends=['ref_profile','merged','scaled_theoretical','deviated']);
    
    '''

rscc_merged = np.array(rscc['interpolated'])
rscc_deviation = np.array(rscc['deviation'])
plt.boxplot(rscc.values(), labels=rscc.keys())
print(rscc_merged.mean())
print(rscc_deviation.mean())