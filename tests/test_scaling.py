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
        from locscale.include.confidenceMapUtil import FDRutil
        self.frequency_map_window = FDRutil.calculate_frequency_map(np.zeros((40, 40, 40)));
        scale_factor_arguments = {}
        scale_factor_arguments['wilson'] = 9.69
        scale_factor_arguments['high_freq'] = 4.69
        scale_factor_arguments['fsc_cutoff'] = 3.4
        scale_factor_arguments['smooth'] = 0.3

        scale_factor_arguments['nyquist'] = 2.6
    
        scale_factor_arguments['boost_secondary_structure'] = 1
        scale_factor_arguments['no_reference'] = False
        self.scale_factor_arguments = scale_factor_arguments
        self.apix = 1.2156
        self.locscale = check_dependencies()['locscale']
        data_folder = os.path.join(self.locscale,'tests','test_data') 
        test_data_dict = {}
        for i in [1,2,3]:
            data_location = os.path.join(data_folder,'test_scaling_data_'+str(i)+'.pickle')
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
            
            print(emmap_window.shape)
            rp_emmap_test,_ = compute_radial_profile(emmap_window, self.frequency_map_window)
            rp_modmap_test, frequencies_map = compute_radial_profile(modmap_window, self.frequency_map_window)
            rp_emmap_target = self.test_data[i]['rp_emmap']
            rp_modmap_target = self.test_data[i]['rp_modmap']
            
            print("Computed radial profiles once")
            
            #self.assertEqual(rp_emmap_test[:-1], rp_emmap_target, msg="Radial profiles of emmap matched target")
            #self.assertEqual(rp_modmap_test[:-1],rp_modmap_target, msg="Radial profiles of modmap matched target")
            
            scale_factors_old,bfactors,qfit = compute_scale_factors(rp_emmap_test, rp_modmap_test, apix=self.apix, scale_factor_arguments=self.scale_factor_arguments, use_theoretical_profile=False)

            scale_factors_new, bfactor,qfit,report = compute_scale_factors(rp_emmap_test, rp_modmap_test, apix=self.apix, scale_factor_arguments=self.scale_factor_arguments, use_theoretical_profile=True, check_scaling=True)
            
            scaled_map_old_test = set_radial_profile(emmap_window, scale_factors_old, frequencies_map, self.frequency_map_window, emmap_window.shape)
            scaled_map_new_test = set_radial_profile(emmap_window, scale_factors_new, frequencies_map, self.frequency_map_window, emmap_window.shape)
            
            scaled_map_old_target = self.test_data[i]['scaled_window_old']
            scaled_map_new_target = self.test_data[i]['scaled_window_new']

    
    def test_scaling_random(self):
        from locscale.utils.scaling_tools import compute_radial_profile, compute_scale_factors, set_radial_profile
        from locscale.include.emmer.ndimage.map_tools import compute_real_space_correlation as rsc
        from locscale.include.emmer.ndimage.map_utils import extract_window
        import mrcfile
        import random
        
        emmap_path = os.path.join(self.locscale,"tests","test_data","emd5778_map_full.mrc")
        modmap_path = os.path.join(self.locscale,"tests","test_data","pdb3j5p_refined_4locscale.mrc")
        print("\nTesting random windows within the maps 100 times")
        for i in range(100):
            emmap = mrcfile.open(emmap_path).data
            modmap = mrcfile.open(modmap_path).data
            
            centers = np.asarray(np.where(emmap>=0.01)).T.tolist()
            random_center_1 = random.choice(centers)
            random_center_2 = random.choice(centers)
            
            emmap_window = extract_window(emmap, random_center_1, 40)
            modmap_window = extract_window(modmap, random_center_2, 40)
            
            rp_emmap_window, _ = compute_radial_profile(emmap_window, self.frequency_map_window)
            
            rp_modmap_window, frequencies_map = compute_radial_profile(modmap_window, self.frequency_map_window)
            
            sf_old,bfactor,qfit = compute_scale_factors(rp_emmap_window, rp_modmap_window, apix=self.apix, 
                                           scale_factor_arguments=self.scale_factor_arguments,
                                           check_scaling=False, use_theoretical_profile=False)
            
            sf_new,bfactor,qfit, report = compute_scale_factors(rp_emmap_window, rp_modmap_window, apix=self.apix, 
                                           scale_factor_arguments=self.scale_factor_arguments,
                                           check_scaling=True, use_theoretical_profile=True)
            
            quality_of_theoretical_profiles = rsc(report['scaled_reference_profile'], report['input_ref_profile'])
            self.assertTrue(quality_of_theoretical_profiles > 0.99)
            scaled_map_old_inside,_ = set_radial_profile(emmap_window, sf_old, frequencies_map, self.frequency_map_window, emmap_window.shape)
            scaled_map_new_inside,_ = set_radial_profile(emmap_window, sf_new, frequencies_map, self.frequency_map_window, emmap_window.shape)
            
            rp_scaled_map_old,_ = compute_radial_profile(scaled_map_old_inside, self.frequency_map_window)
            rp_scaled_map_new,_ = compute_radial_profile(scaled_map_new_inside, self.frequency_map_window)
            
            radial_profile_match_old = rsc(rp_modmap_window,rp_scaled_map_old )
            radial_profile_match_new = rsc(rp_modmap_window,rp_scaled_map_new )
            
            self.assertTrue(radial_profile_match_old > 0.99)
            self.assertTrue(radial_profile_match_new > 0.9)
        
        
        
        
        
        
        
            

            
            
            
            
        
