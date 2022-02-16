#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 10 19:42:23 2021

@author: alok
"""

import unittest
import os

class TestPseudomodelPipeline(unittest.TestCase):
    
    def setUp(self):
        from locscale.pseudomodel.pseudomodel_headers import check_dependencies
        
        self.locscale_path = check_dependencies()['locscale']
        locscale_path = self.locscale_path
        self.emmap_path = os.path.join(locscale_path,"tests","test_data","emd5778_map.mrc")
        self.model_path = os.path.join(locscale_path,"tests","test_data","pdb3j5p_refined.pdb")
        self.mask_path = os.path.join(locscale_path,"tests","test_data","emd5778_mask.mrc")
        self.out_dir = os.path.join(locscale_path,"tests","processed")
        self.wilson_cutoff = 9.69
        self.fsc = 3.4
        self.kick_model = os.path.join(locscale_path,"tests","test_data","pseudomodel.pdb")
    
    def test_pseudomodel_pipeline(self):
        from locscale.pseudomodel.pipeline import get_modmap
        from tempfile import TemporaryDirectory
        
        print("Integration test: Pseudomodel building pipeline")
        
        with TemporaryDirectory() as tempDir: 
            import os
            from subprocess import run
            
            # copy emmap
            run(["cp",self.emmap_path,tempDir])
            temp_emmap_path = os.path.join(tempDir,self.emmap_path.split("/")[-1])
            
            run(["cp",self.mask_path,tempDir])
            temp_mask_path = os.path.join(tempDir, self.mask_path.split("/")[-1])
            
            os.chdir(tempDir)
            modmap_arguments_gradient = {
                'emmap_path':temp_emmap_path,
                'mask_path':temp_mask_path,
                'pdb_path':None,
                'pseudomodel_method':'gradient',
                'pam_distance':1.2,
                'pam_iteration':1,
                'fsc_resolution':3.4,
                'refmac_iter':1,
                'add_blur':0,
                'skip_refine':False,
                'pg_symmetry':"C1",
                'model_resolution':None,
                'molecular_weight':None,
                'build_ca_only':False,
                'verbose':False,
                }
            
            modmap_arguments_kick = {
                'emmap_path':temp_emmap_path,
                'mask_path':temp_mask_path,
                'pdb_path':None,
                'pseudomodel_method':'kick',
                'pam_distance':1.2,
                'pam_iteration':1,
                'fsc_resolution':3.4,
                'refmac_iter':1,
                'add_blur':0,
                'skip_refine':False,
                'pg_symmetry':"C1",
                'model_resolution':None,
                'molecular_weight':None,
                'build_ca_only':False,
                'verbose':False,
                }
            temp_modmap_gradient = get_modmap(modmap_arguments_gradient)
            
            temp_modmap_kick = get_modmap(modmap_arguments_kick)
            
            gradient_modmap_exists = os.path.exists(temp_modmap_gradient)
            kick_modmap_exists = os.path.exists(temp_modmap_kick)
            
            self.assertTrue(gradient_modmap_exists)
            self.assertTrue(kick_modmap_exists)
            
            dir(tempDir)
                   

if __name__ == '__main__':
    unittest.main()           
   
        
        
        
        
    
    
        
        
            
            
        
    
    
        
        
        
        
        
        
        
        
        
        
        
        
