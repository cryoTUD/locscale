#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 13 00:01:47 2021

@author: alok
"""

import unittest
import numpy as np
import os

class test_compute_scaling(unittest.TestCase):
    def setUp(self):
        from locscale.pseudomodel.pseudomodel_headers import check_dependencies
        import pickle
        scale_factor_arguments = {}
        scale_factor_arguments['wilson'] = 8.5
        scale_factor_arguments['high_freq'] = 4.91
        scale_factor_arguments['fsc_cutoff'] = 2.5
        scale_factor_arguments['smooth'] = 0.3
        self.scale_factor_arguments = scale_factor_arguments
        self.apix = 1.2156
        self.locscale = check_dependencies()['locscale']
        data_folder = self.locscale + '/tests/test_data/' 
        test_data_dict = {}
        for i in [1,2,3]:
            data_location = data_folder+'test_scaling_data_'+str(i)+'.pickle'
            with open(data_location,"rb") as f:
                test_data_dict[i] = pickle.load(f)
        
        self.test_data = test_data_dict
        
    
    def test_scaling(self):
        from locscale.utils.scaling_tools import compute_radial_profile, compute_scale_factors, set_radial_profile
        
        for i in [1,2,3]:
            print(i)
            print(self.test_data[i]['report']['scaling_condition'])
            emmap_window = self.test_data[i]['emmap_window']
            modmap_window = self.test_data[i]['modmap_window']
            rp_emmap_test = compute_radial_profile(emmap_window)
            rp_modmap_test, radii = compute_radial_profile(modmap_window, return_indices=True)
            rp_emmap_target = self.test_data[i]['rp_emmap']
            rp_modmap_target = self.test_data[i]['rp_modmap']
            
            self.assertTrue((rp_emmap_test==rp_emmap_target).all())
            self.assertTrue((rp_modmap_test==rp_modmap_target).all())
            
            scale_factors_old = compute_scale_factors(rp_emmap_test, rp_modmap_test, apix=self.apix, scale_factor_arguments=self.scale_factor_arguments, use_theoretical_profile=False)

            scale_factors_new, report = compute_scale_factors(rp_emmap_test, rp_modmap_test, apix=self.apix, scale_factor_arguments=self.scale_factor_arguments, use_theoretical_profile=True, check_scaling=True)
            
            scaled_map_old_test = set_radial_profile(emmap_window, scale_factors_old,radii)
            scaled_map_new_test = set_radial_profile(emmap_window, scale_factors_new, radii)
            
            scaled_map_old_target = self.test_data[i]['scaled_window_old']
            scaled_map_new_target = self.test_data[i]['scaled_window_new']
            print("Not using theoretical profile: ",(scaled_map_old_test==scaled_map_old_target).all())
            print("Using theoretical profile: ",(scaled_map_new_test==scaled_map_new_target).all())
            self.assertTrue((scaled_map_old_test==scaled_map_old_target).all())
            self.assertTrue((scaled_map_new_test==scaled_map_new_target).all())
    
    def test_scaling_random(self):
        from locscale.utils.scaling_tools import compute_radial_profile, compute_scale_factors, set_radial_profile
        from locscale.include.emmer.ndimage.map_tools import compute_real_space_correlation as rsc
        from locscale.include.emmer.ndimage.map_utils import extract_window
        import mrcfile
        import random
        
        emmap_path = self.locscale + "/tests/test_data/emd5778_unfiltered.mrc"
        modmap_path = self.locscale + "/tests/test_data/model_reference.mrc"
        print("Testing random windows within the maps 100 times")
        for i in range(100):
            emmap = mrcfile.open(emmap_path).data
            modmap = mrcfile.open(modmap_path).data
            
            centers = np.asarray(np.where(emmap>=7)).T.tolist()
            random_center_1 = random.choice(centers)
            random_center_2 = random.choice(centers)
            
            emmap_window = extract_window(emmap, random_center_1, 40)
            modmap_window = extract_window(modmap, random_center_2, 40)
            
            rp_emmap_window = compute_radial_profile(emmap_window)
            rp_modmap_window, radii = compute_radial_profile(modmap_window, return_indices=True)
            sf_old = compute_scale_factors(rp_emmap_window, rp_modmap_window, apix=self.apix, 
                                           scale_factor_arguments=self.scale_factor_arguments,
                                           check_scaling=False, use_theoretical_profile=False)
            
            sf_new, report = compute_scale_factors(rp_emmap_window, rp_modmap_window, apix=self.apix, 
                                           scale_factor_arguments=self.scale_factor_arguments,
                                           check_scaling=True, use_theoretical_profile=True)
            
            quality_of_theoretical_profiles = rsc(report['scaled_reference_profile'], report['input_ref_profile'])
            self.assertTrue(quality_of_theoretical_profiles > 0.99)
            
            rp_scaled_map_old = compute_radial_profile(set_radial_profile(emmap_window, sf_old, radii))
            rp_scaled_map_new = compute_radial_profile(set_radial_profile(emmap_window, sf_new, radii))
            
            radial_profile_match_old = rsc(rp_modmap_window,rp_scaled_map_old )
            radial_profile_match_new = rsc(rp_modmap_window,rp_scaled_map_new )
            
            self.assertTrue(radial_profile_match_old > 0.99)
            self.assertTrue(radial_profile_match_new > 0.99)
        
        
        
        
        
        
        
            

            
            
            
            
        