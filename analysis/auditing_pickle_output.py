#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Sep 11 11:17:21 2021

@author: alok
"""
import pickle
import numpy as np
import matplotlib.pyplot as plt
from emmer.ndimage.profile_tools import plot_radial_profile

pickle_output = "/mnt/c/Users/abharadwaj1/Downloads/ForUbuntu/LocScale/emd5778/dst_1o2_using_zsh/profiles_audit_after_proper_resample.pickle"

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
with open(pickle_output,"rb") as audit_file:
    audit_scaling = pickle.load(audit_file)

random_positions = list(audit_scaling.keys())    
for key in random_positions[:1]:
    freq = audit_scaling[key]['freq']
    em_profile = audit_scaling[key]['em_profile']
    ref_profile = audit_scaling[key]['input_ref_profile']
    theoretical_profile = audit_scaling[key]['theoretical_amplitude']
    scaled_theoretical = audit_scaling[key]['scaled_theoretical_amplitude']
    merged_profile = audit_scaling[key]['scaled_reference_profile']
    
    plot_radial_profile(freq,[em_profile, ref_profile, theoretical_profile, scaled_theoretical, merged_profile],legends=['em_profile','ref_profile','th profile','scaled th profile','merged'])
    
    