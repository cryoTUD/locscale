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
        import os
        self.locscale_path = get_locscale_path()
        
        self.emmap_path = os.path.join(self.locscale_path,"tests","test_data","emd5778_map_full.mrc")
                
        
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
   
        
        
        
        
    
    
        
        
            
            
        
    
    
        
        
        
        
        
        
        
        
        
        
        
        
