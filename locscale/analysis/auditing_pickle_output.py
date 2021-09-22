#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Sep 11 11:17:21 2021

@author: alok
"""
import pickle
import numpy as np
import matplotlib.pyplot as plt
from locscale.include.emmer.ndimage.profile_tools import plot_radial_profile

folder = "/mnt/c/Users/abharadwaj1/Downloads/ForUbuntu/LocScale/emd5778/"

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

pickle_output = folder+"profiles_audit.pickle"
with open(pickle_output,"rb") as audit_file:
    audit_scaling = pickle.load(audit_file)

import random
random_positions = list(audit_scaling.keys())    
sample = random.sample(random_positions, 1)
for key in sample:
    freq = audit_scaling[key]['freq']
    em_profile = audit_scaling[key]['em_profile']
    ref_profile = audit_scaling[key]['input_ref_profile']
    theoretical_profile = audit_scaling[key]['theoretical_amplitude']
    scaled_theoretical = audit_scaling[key]['scaled_theoretical_amplitude']
    merged_profile = audit_scaling[key]['scaled_reference_profile']
    print(audit_scaling[key]['scaling_condition'])
    
    
    plot_radial_profile(freq,[em_profile, ref_profile, theoretical_profile, scaled_theoretical, merged_profile],legends=['em_profile','ref_profile','th profile','scaled th profile','merged']);
    
    