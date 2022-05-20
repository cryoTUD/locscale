#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 13 00:01:47 2021

@author: alok
"""

import unittest
import numpy as np
import os

class test_locscale(unittest.TestCase):
    def setUp(self):
        from locscale.utils.file_tools import get_locscale_path
        import pickle
        from locscale.include.confidenceMapUtil import FDRutil

        self.locscale = get_locscale_path()
        data_folder = os.path.join(self.locscale,'tests','test_data') 
        self.emmap_path = os.path.join(data_folder, "emd5778_map_chainA.mrc")
        self.mask_path = os.path.join(data_folder, "emd5778_mask_chainA.mrc")
        self.model_coordinates = os.path.join(data_folder, "pdb3j5p_refined_chainA.pdb")
        self.reference_locscale_MB = os.path.join(data_folder, "reference_mb_locscale.mrc")
        self.reference_locscale_MF = os.path.join(data_folder, "reference_mf_locscale.mrc")
        self.resolution = 3.4
        
        
    
    def test_run_model_based_locscale(self):
        from tempfile import TemporaryDirectory
        
        print("Testing: Model Based LocScale")
        with TemporaryDirectory() as tempDir: 
            from locscale.include.emmer.ndimage.map_tools import compute_real_space_correlation as rsc
            import os
            from subprocess import run
            
            # copy emmap
            run(["cp",self.emmap_path,tempDir])
            emmap_name= os.path.basename(self.emmap_path)
            temp_emmap_path = os.path.join(tempDir,emmap_name)
            
            run(["cp",self.model_coordinates,tempDir])
            model_name = os.path.basename(self.model_coordinates)
            temp_model_path = os.path.join(tempDir, model_name)

            run(["cp",self.mask_path,tempDir])
            mask_name = os.path.basename(self.mask_path)
            temp_mask_path = os.path.join(tempDir, mask_name)
            
            os.chdir(tempDir)
            
            output_locscale_path = os.path.join(tempDir, "locscale_unittest.mrc")
            locscale_script_path = os.path.join(self.locscale,"locscale","main.py")
            
            locscale_command = ["python",locscale_script_path,"run_locscale","--emmap_path",\
                temp_emmap_path, "--model_coordinates",temp_model_path,"--mask",temp_mask_path, \
                "--ref_resolution","3.4","--outfile",output_locscale_path,"--skip_refine","--verbose"]
            
            locscale_test_run = run(locscale_command)
            
            self.assertTrue(os.path.exists(output_locscale_path))
            
            rscc_test = rsc(self.reference_locscale_MB,output_locscale_path)
            
            self.assertTrue(rscc_test>0.99)
            
    def test_run_model_free_locscale(self):
        from tempfile import TemporaryDirectory
        
        print("Testing: Model Free LocScale")
        with TemporaryDirectory() as tempDir: 
            from locscale.include.emmer.ndimage.map_tools import compute_real_space_correlation as rsc
            import os
            from subprocess import run
            
            # copy emmap
            run(["cp",self.emmap_path,tempDir])
            emmap_name= os.path.basename(self.emmap_path)
            temp_emmap_path = os.path.join(tempDir,emmap_name)

            run(["cp",self.mask_path,tempDir])
            mask_name = os.path.basename(self.mask_path)
            temp_mask_path = os.path.join(tempDir, mask_name)
            
            os.chdir(tempDir)
            
            output_locscale_path = os.path.join(tempDir, "locscale_MF_unittest.mrc")
            locscale_script_path = os.path.join(self.locscale,"locscale","main.py")
            
            locscale_command = ["python",locscale_script_path,"run_locscale","--emmap_path",temp_emmap_path, \
                "--mask",temp_mask_path, "--outfile",output_locscale_path,"--ref_resolution","3.4","--verbose"]
                        
            locscale_test_run = run(locscale_command)
            
            self.assertTrue(os.path.exists(output_locscale_path))
            
            rscc_test = rsc(self.reference_locscale_MF,output_locscale_path)
            
            self.assertTrue(rscc_test>0.99)
            
    
   
            
            
            
            
        
