#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 10 19:42:23 2021

@author: alok
"""

import unittest

      

class TestSymmetry(unittest.TestCase):
    
    def setUp(self):
        from locscale.utils.file_tools import get_locscale_path
        
        self.locscale_path = get_locscale_path()
        lPath = self.locscale_path
        self.emmap_path = lPath+"/tests/test_data/emd5778_map.mrc"
        self.model_path = lPath+"/tests/test_data/pdb3j5p_refined.pdb"
        self.mask_path = lPath+"/tests/test_data/emd5778_mask.mrc"
        self.out_dir = lPath+"/tests/processed/"
        self.wilson_cutoff = 8.55
        self.fsc = 3.4
        self.kick_model = lPath+"/tests/test_data/pseudomodel.pdb"
        
                
        
    def test_symmetry(self):
        print("Imposing a symmetry condition of C4")
        import emda.emda_methods as em
        from tempfile import TemporaryDirectory
        
        with TemporaryDirectory() as tempDir:
            import os
            os.chdir(tempDir)
            sym = em.symmetry_average([self.emmap_path],[3.4],pglist=["C4"])
            self.assertEqual(sym[0].shape,(256,256,256))
        
        
    
    
        

if __name__ == '__main__':
    unittest.main()           
   
        
        
        
        
    
    
        
        
            
            
        
    
    
        
        
        
        
        
        
        
        
        
        
        
        
