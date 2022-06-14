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
        
        self.emmap_path = os.path.join(self.locscale_path,"locscale","tests","test_data","emd5778_map_chainA.mrc")
        #self.emmap_path = os.path.join(self.locscale_path,"tests","test_data","emd5778_map_full.mrc")
        self.symmetry_output = os.path.join(self.locscale_path,"locscale","tests","test_data","avgmap_fcodes.mrc")
#        self.emmap_path = "/home/alok/dev/ForUbuntu/LocScale/tests/new_symmetry/emd5778_tutorial.mrc"

    def copy_files(self, file_path, tempDir):
        import os
        from subprocess import run
        run(["cp",file_path,tempDir])
        new_file_path = os.path.join(tempDir, os.path.basename(file_path))
        if os.path.exists(new_file_path):        
            from locscale.include.emmer.ndimage.map_utils import load_map, save_as_mrc, resample_map
            emmap, apix = load_map(new_file_path)   
            resampled_emmap=resample_map(emmap, apix=apix, apix_new=3)
            resampled_file_path = os.path.join(tempDir, "resampled.mrc")
            print("Saving resampled map to: {}".format(resampled_file_path))
            print("Resampled map shape: {}".format(resampled_emmap.shape))
            save_as_mrc(resampled_emmap, resampled_file_path, apix=3) 
            return resampled_file_path
            
        else:
            raise UserWarning("Could not copy {} to {}".format(file_path,tempDir))

    def test_symmetry(self):
        print("Imposing a symmetry condition of C4")
        from locscale.include.symmetry_emda.symmetrize_map import symmetrize_map_emda
        from tempfile import TemporaryDirectory
        import os
        from locscale.include.emmer.ndimage.map_utils import load_map, save_as_mrc, resample_map
        from locscale.include.emmer.ndimage.map_tools import compute_real_space_correlation as rsc

        with TemporaryDirectory() as tempDir:
            copied_emmap_path = self.copy_files(self.emmap_path, tempDir)
            os.chdir(tempDir)
            emmap, apix = load_map(copied_emmap_path)  
            print("Emmap map shape: {}".format(emmap.shape))         
            sym = symmetrize_map_emda(copied_emmap_path, pg="C4")
         #   self.assertEqual(sym.shape,(104,104,104))
            reference_symmetry_map, _ = load_map(self.symmetry_output)
            rscc = rsc(sym, reference_symmetry_map)
            print("Correlation with test map", rscc)
            self.assertTrue(rscc>0.99)
            
            
if __name__ == '__main__':
    unittest.main()           
   
        
        
        
        
    
    
        
        
            
            
        
    
    
        
        
        
        
        
        
        
        
        
        
        
        
