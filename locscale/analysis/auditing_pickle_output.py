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

folder = "/mnt/c/Users/abharadwaj1/Downloads/ForUbuntu/LocScale/tests/map_sharpening/Milan_APO/using_deviations_addition"

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

pickle_output = os.path.join(folder,"profiles_audit.pickle")

import mrcfile
#model_map = mrcfile.open("/mnt/c/Users/abharadwaj1/Downloads/ForUbuntu/LocScale/tests/test_low_frequency_dip/pdb3j5p_refined_cropped_shifted_refmac_refined_4locscale.mrc").data

with open(pickle_output,"rb") as audit_file:
    audit_scaling = pickle.load(audit_file)

import random
random_positions = list(audit_scaling.keys())    
sample = random.sample(random_positions, 1)
rscc = {
        'interpolated':[],
        'deviation':[]}

deviation_magnitudes = {}
for key in sample:
#    model_map_window = get_box(model_map, key, size=24)
    
    freq = audit_scaling[key]['freq']
 #   rp_model_map = compute_radial_profile(model_map_window)
    em_profile = audit_scaling[key]['em_profile']
    ref_profile = audit_scaling[key]['input_ref_profile']
    theoretical_profile = audit_scaling[key]['theoretical_amplitude']
    scaled_theoretical = audit_scaling[key]['scaled_theoretical_amplitude']
    merged_profile = audit_scaling[key]['scaled_reference_profile']
    deviated_profile = audit_scaling[key]['deviated_reference_profile']


    bfactor = (np.log(audit_scaling[key]['exponential_fit'][0])-np.log(audit_scaling[key]['exponential_fit'][-1])) / (freq[0]**2 - freq[-1]**2)
    deviation = deviated_profile / ref_profile
    
    deviation_magnitudes[key] = [bfactor, deviation.mean()]
    '''
    print(audit_scaling[key]['scaling_condition'])
    
    plot_radial_profile(freq,[ref_profile],
                        legends=['ref_profile']);
    
    plot_radial_profile(freq,[scaled_theoretical],
                        legends=['scaled theoretical']);
    
    plot_radial_profile(freq,[ref_profile, scaled_theoretical],
                        legends=['ref_profile','scaled theoretical']);
    '''    
    plot_radial_profile(freq,[ref_profile, merged_profile, scaled_theoretical, deviated_profile, merged_profile],
                        legends=['ref_profile','merged','scaled_theoretical','deviated', 'merged_profile']);
    


