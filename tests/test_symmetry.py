#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 10 19:42:23 2021

@author: alok
"""

import unittest
import numpy as np
import os

      

class TestPseudomodelHeaders(unittest.TestCase):
    
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
        
                
        
    def test_symmetry(self):
        print("Imposing a symmetry condition of C4")
        import locscale.include.emda.emda.emda_methods as em
        from tempfile import TemporaryDirectory
        
        with TemporaryDirectory() as tempDir:
            sym = em.symmetry_average([self.emmap_path],[3.4],pglist=["C4"])
            self.assertEqual(sym[0].shape,(256,256,256))
        
        
    
    
        

if __name__ == '__main__':
    unittest.main()           
   
        
        
        
        
    
    
        
        
            
            
        
    
    
        
        
        
        
        
        
        
        
        
        
        
        