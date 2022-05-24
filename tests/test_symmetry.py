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

    def copy_files(self, file_path, tempDir):
        from subprocess import run
        run(["cp",file_path,tempDir])
        if os.path.exists(os.path.join(tempDir, os.path.basename(file_path))):               
            return os.path.join(tempDir, os.path.basename(file_path))
            
        else:
            raise UserWarning("Could not copy {} to {}".format(path,tempDir))

        
    def test_symmetry(self):
        print("Imposing a symmetry condition of C4")
        import emda.emda_methods as em
        from tempfile import TemporaryDirectory
        
        with TemporaryDirectory() as tempDir:
            copied_emmap_path = self.copy_files(self.emmap_path, tempDir)
            sym = em.symmetry_average([copied_emmap_path],[3.4],pglist=["C4"])
            self.assertEqual(sym[0].shape,(256,256,256))
        
       

if __name__ == '__main__':
    unittest.main()           
   
        
        
        
        
    
    
        
        
            
            
        
    
    
        
        
        
        
        
        
        
        
        
        
        
        
