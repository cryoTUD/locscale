#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 10 19:42:23 2021

@author: alok
"""

import unittest
     

class TestPseudomodelPipeline(unittest.TestCase):
    
    def setUp(self):
        from locscale.pseudomodel.pseudomodel_headers import check_dependencies
        
        self.locscale_path = check_dependencies()['locscale']
        lPath = self.locscale_path
        self.emmap_path = lPath+"/tests/test_data/emd5778_unfiltered.mrc"
        self.model_path = lPath+"/tests/test_data/pdb3j5p_refined_cropped.pdb"
        self.mask_path = lPath+"/tests/test_data/emd5778_mask.mrc"
        self.out_dir = lPath+"/tests/processed/"
        self.wilson_cutoff = 8.55
        self.fsc = 3.4
        self.kick_model = lPath+"/tests/test_data/kick_pseudomodel.pdb"
    
    def test_pseudomodel_pipeline(self):
        from locscale.pseudomodel.pipeline import get_modmap
        from tempfile import TemporaryDirectory
        
        print("Integration test: Pseudomodel building pipeline")
        
        with TemporaryDirectory() as tempDir: 
            import os
            from subprocess import run
            
            # copy emmap
            run(["cp",self.emmap_path,tempDir])
            temp_emmap_path = tempDir + "/" + self.emmap_path.split("/")[-1]
            
            run(["cp",self.mask_path,tempDir])
            temp_mask_path = tempDir + "/" + self.mask_path.split("/")[-1]
            
            os.chdir(tempDir)
            temp_modmap_gradient = get_modmap(
                temp_emmap_path, temp_mask_path, pdb_path=None,
                pseudomodel_method='gradient',pam_distance=1.2, 
                pam_iteration=1, fsc_resolution=3.4, 
                refmac_iter=1, add_blur=0, verbose=True)
            
            temp_modmap_kick = get_modmap(
                temp_emmap_path, temp_mask_path, pdb_path=None,
                pseudomodel_method='kick',pam_distance=1.2, 
                pam_iteration=1, fsc_resolution=3.4, 
                refmac_iter=1, add_blur=0, verbose=True)
            
            gradient_modmap_exists = os.path.exists(temp_modmap_gradient)
            kick_modmap_exists = os.path.exists(temp_modmap_kick)
            
            self.assertTrue(gradient_modmap_exists)
            self.assertTrue(kick_modmap_exists)
            
            dir(tempDir)
            
            
            
        

    
    
        

if __name__ == '__main__':
    unittest.main()           
   
        
        
        
        
    
    
        
        
            
            
        
    
    
        
        
        
        
        
        
        
        
        
        
        
        