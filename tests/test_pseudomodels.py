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
        from scripts.get_pseudomodel.pseudomodel_headers import check_dependencies
        
        self.locscale_path = check_dependencies()['locscale']
        lPath = self.locscale_path
        self.emmap_path = lPath+"/tests/test_data/emd5778_unfiltered.mrc"
        self.model_path = lPath+"/tests/test_data/pdb3j5p_refined_cropped.pdb"
        self.mask_path = lPath+"/tests/test_data/emd5778_mask.mrc"
        self.out_dir = lPath+"/tests/processed/"
        self.wilson_cutoff = 8.55
        self.fsc = 3.4
        self.kick_model = lPath+"/tests/test_data/kick_pseudomodel.pdb"
        
                
        
    def test_sharpen_maps(self):
        from scripts.get_pseudomodel.pseudomodel_headers import prepare_sharpen_map
        
        print("Testing: prepare_sharpen_map \n")
        outfile, pwlf_fit = prepare_sharpen_map(emmap_path=self.emmap_path, wilson_cutoff=self.wilson_cutoff, fsc_resolution=self.fsc, return_processed_files=True)
        
        sharpened_map_present = os.path.exists(outfile)
        self.assertTrue(sharpened_map_present)
        
        slopes = pwlf_fit.calc_slopes()
        self.assertEqual(len(slopes),3)
        self.assertTrue(slopes[0]<0 and slopes[1] > 0 and slopes[2] < 0)
        self.assertAlmostEqual(slopes[2]*4, -61, delta=5.0)
        
        f2_breakpoints = pwlf_fit.fit(n_segments=3)
        d_breakpoints = np.sqrt(1/f2_breakpoints)
        self.assertTrue(d_breakpoints[0] < 8.55 and d_breakpoints[1] > 5 and d_breakpoints[2] > 4)
        
        ## Remove files
        os.remove(outfile)
        
    def test_run_FDR(self):
        from scripts.get_pseudomodel.pseudomodel_headers import run_FDR
        print("Testing: run_FDR")
        import mrcfile
        mask_path = run_FDR(emmap_path=self.emmap_path, window_size=40, verbose=False)
        
        mask_exists = os.path.exists(mask_path)
        self.assertTrue(mask_exists)
        mask = mrcfile.open(mask_path).data
        self.assertAlmostEqual(mask.sum(), 277867, delta=27786.7)
        
        os.remove(mask_path)
        
    def test_run_pam(self):
        from scripts.get_pseudomodel.pseudomodel_headers import run_pam
        import gemmi
        
        def quick_check_pseudomodel(pseudomodel_path):
            gemmi_st = gemmi.read_structure(pseudomodel_path_gradient)
        
            num_atoms_final = gemmi_st[0].count_atom_sites()
            center_of_mass = np.array(gemmi_st[0].calculate_center_of_mass().tolist())
            displacement = np.array([256*1.2156/2,256*1.2156/2,256*1.2156/2]) - center_of_mass
            distance_from_center = np.linalg.norm(displacement)
            return num_atoms_final, distance_from_center
        
        print("Testing: run_pam for gradient method")
        pseudomodel_path_gradient = run_pam(emmap_path=self.emmap_path,mask_path=self.mask_path,threshold=1,num_atoms=16040,method='gradient',bl=1.2,g=None,friction=None,scale_map=None,scale_lj=None,total_iterations=3,verbose=False)
        
        gradient_pseudomodel_exists = os.path.exists(pseudomodel_path_gradient)
        self.assertTrue(gradient_pseudomodel_exists)
        
        num_atoms_g, center_offset_g = quick_check_pseudomodel(pseudomodel_path_gradient)
        
        self.assertEqual(num_atoms_g, 16040)
        self.assertLess(center_offset_g, 50)
        
        print("Testing: run_pam with kick method")
        pseudomodel_path_kick = run_pam(emmap_path=self.emmap_path,mask_path=self.mask_path,threshold=1,num_atoms=16040,method='kick',bl=1.2,g=None,friction=None,scale_map=None,scale_lj=None,total_iterations=3,verbose=False)
        
        kick_pseudomodel_exists = os.path.exists(pseudomodel_path_kick)
        self.assertTrue(kick_pseudomodel_exists)
        
        num_atoms_k, center_offset_k = quick_check_pseudomodel(pseudomodel_path_kick)
        
        self.assertEqual(num_atoms_k, 16040)
        self.assertLess(center_offset_k, 50)
        
        os.remove(pseudomodel_path_kick)
        os.remove(pseudomodel_path_gradient)
        
    def test_run_refmac(self):
        from scripts.get_pseudomodel.pseudomodel_headers import run_refmac
        print("Testing: run_refmac refinement")
        refined_model=run_refmac(model_path=self.kick_model,model_name=self.kick_model[:-4],map_path=self.emmap_path,resolution=self.fsc,maskdims=[256*1.2156,256*1.2156,256*1.2156],  num_iter=1,verbose=False)
        
        refined_model_path_exists = os.path.exists(refined_model)
        self.assertTrue(refined_model_path_exists)
        
        os.remove(refined_model)
    
    def test_run_refmap(self):
        from scripts.get_pseudomodel.pseudomodel_headers import run_refmap
        print("Testing: run_refmap")
        refmap_model = run_refmap(model_path=self.model_path,emmap_path=self.emmap_path,mask_path=self.mask_path,resolution=self.fsc,verbose=False)
        
        refmap_pseudomodel = run_refmap(model_path=self.kick_model,emmap_path=self.emmap_path,mask_path=self.mask_path,resolution=self.fsc,verbose=False)
        
        self.assertTrue(os.path.exists(refmap_model))
        self.assertTrue(os.path.exists(refmap_pseudomodel))
        
        print("Testing if mapmask.sh is present")
        
        mapmask_present = os.path.exists(self.locscale_path + "/scripts/utils/mapmask.sh")
        self.assertTrue(mapmask_present)
        
        os.remove(refmap_pseudomodel)
        os.remove(refmap_model)
        

if __name__ == '__main__':
    unittest.main()           
   
        
        
        
        
    
    
        
        
            
            
        
    
    
        
        
        
        
        
        
        
        
        
        
        
        