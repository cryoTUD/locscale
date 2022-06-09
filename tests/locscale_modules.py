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
        from locscale.utils.file_tools import get_locscale_path
        import os
        
        self.locscale_path = get_locscale_path()
        
        self.emmap_path_full = os.path.join(self.locscale_path,"tests","test_data","emd5778_map_full.mrc")
        self.model_path_full = os.path.join(self.locscale_path,"tests","test_data","pdb3j5p_refined.pdb")
        
        self.mask_path_full = os.path.join(self.locscale_path,"tests","test_data","emd5778_mask_full.mrc")
        self.out_dir = os.path.join(self.locscale_path,"tests","processed")
        self.wilson_cutoff = 9.69
        self.fsc = 3.4
        self.pseudomodel_full = os.path.join(self.locscale_path,"tests","test_data","pseudomodel.pdb")
        
    
    def copy_files(self, path, tempDir):
            from subprocess import run
            run(["cp",path,tempDir])
            if os.path.exists(os.path.join(tempDir, os.path.basename(path))):
                return os.path.join(tempDir, os.path.basename(path))
            
            else:
                raise UserWarning("Could not copy {} to {}".format(path,tempDir))
        
                
        
    def test_sharpen_maps(self):
        from locscale.preprocessing.headers import prepare_sharpen_map
        from tempfile import TemporaryDirectory
        print("Testing: prepare_sharpen_map \n")
        
        
        with TemporaryDirectory() as tempDir: 
            from subprocess import run
            temp_emmap_path = self.copy_files(self.emmap_path_full,tempDir)
            os.chdir(tempDir)
            
            outfile, pwlf_fit = prepare_sharpen_map(emmap_path=temp_emmap_path, wilson_cutoff=self.wilson_cutoff, fsc_resolution=self.fsc, return_processed_files=True)
            
            sharpened_map_present = os.path.exists(outfile)
            self.assertTrue(sharpened_map_present)
            
            slopes = pwlf_fit.calc_slopes()
            self.assertEqual(len(slopes),3)
            self.assertTrue(slopes[0]<0 and slopes[1] > 0 and slopes[2] < 0)
            self.assertAlmostEqual(slopes[2]*4, -221, delta=5.0)
            
            f2_breakpoints = pwlf_fit.fit(n_segments=3)
            d_breakpoints = np.sqrt(1/f2_breakpoints)
            self.assertAlmostEqual(d_breakpoints[0], 9.5, delta=0.1)
            self.assertAlmostEqual(d_breakpoints[1], 6.2, delta=0.1)
            self.assertAlmostEqual(d_breakpoints[2], 4.7, delta=0.1)
            self.assertAlmostEqual(d_breakpoints[3], 3.4, delta=0.1)
            
        
    def test_run_FDR(self):
        from locscale.preprocessing.headers import run_FDR
        from tempfile import TemporaryDirectory
        print("Testing: run_FDR")
        import mrcfile
        
        with TemporaryDirectory() as tempDir: 
            temp_emmap_path = self.copy_files(self.emmap_path_full, tempDir)
            os.chdir(tempDir)
            mask_path = run_FDR(emmap_path=temp_emmap_path, window_size=40, verbose=False)
            mask_exists = os.path.exists(mask_path)
            self.assertTrue(mask_exists)
            mask = mrcfile.open(mask_path).data
            self.assertAlmostEqual(mask.sum(), 414495, delta=200)
            
        
    def test_run_pam(self):
        from locscale.preprocessing.headers import run_pam
        import gemmi
        from tempfile import TemporaryDirectory
    
            
        with TemporaryDirectory() as tempDir: 
            temp_emmap = self.copy_files(self.emmap_path_full, tempDir)
            temp_mask = self.copy_files(self.mask_path_full, tempDir)
            os.chdir(tempDir)
                       
        
            def quick_check_pseudomodel(pseudomodel_path):
                gemmi_st = gemmi.read_structure(pseudomodel_path_gradient)
            
                num_atoms_final = gemmi_st[0].count_atom_sites()
                center_of_mass = np.array(gemmi_st[0].calculate_center_of_mass().tolist())
                displacement = np.array([256*1.2156/2,256*1.2156/2,256*1.2156/2]) - center_of_mass
                distance_from_center = np.linalg.norm(displacement)
                return num_atoms_final, distance_from_center
            
            print("Testing: run_pam for gradient method")
            pseudomodel_path_gradient = run_pam(emmap_path=temp_emmap,mask_path=temp_mask,threshold=1,num_atoms=16040,method='gradient',bl=1.2,g=None,friction=None,scale_map=None,scale_lj=None,total_iterations=3,verbose=False)
            
            gradient_pseudomodel_exists = os.path.exists(pseudomodel_path_gradient)
            self.assertTrue(gradient_pseudomodel_exists)
            
            num_atoms_g, center_offset_g = quick_check_pseudomodel(pseudomodel_path_gradient)
            
            self.assertEqual(num_atoms_g, 16040)
            self.assertLess(center_offset_g, 50)
            
            print("Testing: run_pam with kick method")
            pseudomodel_path_kick = run_pam(emmap_path=temp_emmap,mask_path=temp_mask,threshold=1,num_atoms=16040,method='kick',bl=1.2,g=None,friction=None,scale_map=None,scale_lj=None,total_iterations=3,verbose=False)
            
            kick_pseudomodel_exists = os.path.exists(pseudomodel_path_kick)
            self.assertTrue(kick_pseudomodel_exists)
            
            num_atoms_k, center_offset_k = quick_check_pseudomodel(pseudomodel_path_kick)
            
            self.assertEqual(num_atoms_k, 16040)
            self.assertLess(center_offset_k, 50)
            
           
        
    def test_run_servalcat(self):
        from locscale.preprocessing.headers import run_refmac_servalcat
        from tempfile import TemporaryDirectory
        
        print("Testing: run_refmac refinement")
        with TemporaryDirectory() as tempDir: 
            import os
            from subprocess import run
            
            temp_emmap_path = self.copy_files(self.emmap_path_full, tempDir)
            temp_model_path = self.copy_files(self.model_path_full, tempDir)
            
            
            os.chdir(tempDir)
            refined_model=run_refmac_servalcat(model_path=temp_model_path,map_path=temp_emmap_path,resolution=self.fsc,num_iter=1,only_bfactor_refinement=True, verbose=True)
            
            refined_model_path_exists = os.path.exists(refined_model)
            self.assertTrue(refined_model_path_exists)
        
        
    
    def test_run_refmap(self):
        from locscale.preprocessing.headers import run_refmap
        from tempfile import TemporaryDirectory
        print("Testing: run_refmap")
        with TemporaryDirectory() as tempDir: 
            import os
            from subprocess import run
            temp_model_path = self.copy_files(self.model_path_full, tempDir)
            temp_emmap_path = self.copy_files(self.emmap_path_full, tempDir)
            temp_pseudomodel_path = self.copy_files(self.pseudomodel_full, tempDir)
            temp_mask_path = self.copy_files(self.mask_path_full, tempDir)
            
            os.chdir(tempDir)
            refmap_model = run_refmap(model_path=temp_model_path,emmap_path=temp_emmap_path,mask_path=temp_mask_path,resolution=self.fsc,verbose=False)
        
            refmap_pseudomodel = run_refmap(model_path=temp_pseudomodel_path,emmap_path=temp_emmap_path,mask_path=temp_mask_path,resolution=self.fsc,verbose=False)
        
            self.assertTrue(os.path.exists(refmap_model))
            self.assertTrue(os.path.exists(refmap_pseudomodel))
            
            print("Testing if mapmask.sh is present")
            
            mapmask_present = os.path.exists(self.locscale_path + "/locscale/utils/mapmask.sh")
            self.assertTrue(mapmask_present)
            

        
    
        

if __name__ == '__main__':
    unittest.main()           
   
        
        
        
        
    
    
        
        
            
            
        
    
    
        
        
        
        
        
        
        
        
        
        
        
        