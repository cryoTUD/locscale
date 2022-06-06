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
        self.emmap_path = "/home/alok/dev/ForUbuntu/LocScale/tests/new_symmetry/emd5778_tutorial.mrc"

    def copy_files(self, file_path, tempDir):
        import os
        from subprocess import run
        run(["cp",file_path,tempDir])
        if os.path.exists(os.path.join(tempDir, os.path.basename(file_path))):               
            return os.path.join(tempDir, os.path.basename(file_path))
            
        else:
            raise UserWarning("Could not copy {} to {}".format(file_path,tempDir))

        
    def test_symmetry(self):
        print("Imposing a symmetry condition of C4")
        from locscale.include.symmetry_emda.symmetrize_map import symmetrize_map_known_pg
        from tempfile import TemporaryDirectory
        import os
        from locscale.include.emmer.ndimage.map_utils import load_map, save_as_mrc, resample_map

        with TemporaryDirectory() as tempDir:
            copied_emmap_path = self.copy_files(self.emmap_path, tempDir)
            emmap, apix = load_map(copied_emmap_path)
            resampled_emmap = resample_map(emmap, apix=apix,apix_new=3)
            resampled_emmap_path = os.path.join(tempDir, "resampled_emmap.mrc")
            save_as_mrc(resampled_emmap, resampled_emmap_path, apix=3)
            os.chdir(tempDir)
            sym = symmetrize_map_known_pg(resampled_emmap, apix=3, pg="C4")
            self.assertEqual(sym.shape,(104,104,104))
            
        
       

if __name__ == '__main__':
    unittest.main()           
   
        
        
        
        
    
    
        
        
            
            
        
    
    
        
        
        
        
        
        
        
        
        
        
        
        
