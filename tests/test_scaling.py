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
        from pseudomodel_headers import check_dependencies
        import pickle
        scale_factor_arguments = {}
        scale_factor_arguments['wilson'] = 8.55
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
        from compute import compute_radial_profile, compute_scale_factors, set_radial_profile
        
        test_file = {}
        for i in [1,2,3]:
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
            
            self.assertTrue((scaled_map_old_test==scaled_map_old_target).all())
            self.assertTrue((scaled_map_new_test==scaled_map_new_target).all())
            
            
            
        