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
        self.emmap_path = lPath+"/tests/test_data/emd5778_map.mrc"
        self.model_path = lPath+"/tests/test_data/pdb3j5p_refined.pdb"
        self.mask_path = lPath+"/tests/test_data/emd5778_mask.mrc"
        self.out_dir = lPath+"/tests/processed/"
        self.wilson_cutoff = 9.69
        self.fsc = 3.4
        self.kick_model = lPath+"/tests/test_data/pseudomodel.pdb"
        
                
        
    def test_sharpen_maps(self):
        from locscale.pseudomodel.pseudomodel_headers import prepare_sharpen_map
        
        print("Testing: prepare_sharpen_map \n")
        outfile, pwlf_fit = prepare_sharpen_map(emmap_path=self.emmap_path, wilson_cutoff=self.wilson_cutoff, fsc_resolution=self.fsc, return_processed_files=True)
        
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
        
        ## Remove files
        os.remove(outfile)
        
    def test_run_FDR(self):
        from locscale.pseudomodel.pseudomodel_headers import run_FDR
        print("Testing: run_FDR")
        import mrcfile
        mask_path = run_FDR(emmap_path=self.emmap_path, window_size=40, verbose=False)
        
        mask_exists = os.path.exists(mask_path)
        self.assertTrue(mask_exists)
        mask = mrcfile.open(mask_path).data
        self.assertAlmostEqual(mask.sum(), 414495, delta=200)
        
        os.remove(mask_path)
        
    def test_run_pam(self):
        from locscale.pseudomodel.pseudomodel_headers import run_pam
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
        from locscale.pseudomodel.pseudomodel_headers import run_refmac
        from tempfile import TemporaryDirectory
        
        print("Testing: run_refmac refinement")
        with TemporaryDirectory() as tempDir: 
            import os
            from subprocess import run
            
            # copy emmap
            run(["cp",self.emmap_path,tempDir])
            temp_emmap_path = tempDir + "/" + self.emmap_path.split("/")[-1]
            
            run(["cp",self.model_path,tempDir])
            temp_model_path = tempDir + "/" + self.model_path.split("/")[-1]
            os.chdir(tempDir)
            refined_model=run_refmac(model_path=temp_model_path,map_path=temp_emmap_path,resolution=self.fsc,num_iter=1,only_bfactor_refinement=True, verbose=True)
            
            refined_model_path_exists = os.path.exists(refined_model)
            self.assertTrue(refined_model_path_exists)
        
        
    
    def test_run_refmap(self):
        from locscale.pseudomodel.pseudomodel_headers import run_refmap
        print("Testing: run_refmap")
        refmap_model = run_refmap(model_path=self.model_path,emmap_path=self.emmap_path,mask_path=self.mask_path,resolution=self.fsc,verbose=False)
        
        refmap_pseudomodel = run_refmap(model_path=self.kick_model,emmap_path=self.emmap_path,mask_path=self.mask_path,resolution=self.fsc,verbose=False)
        
        self.assertTrue(os.path.exists(refmap_model))
        self.assertTrue(os.path.exists(refmap_pseudomodel))
        
        print("Testing if mapmask.sh is present")
        
        mapmask_present = os.path.exists(self.locscale_path + "/locscale/utils/mapmask.sh")
        self.assertTrue(mapmask_present)
        
        os.remove(refmap_pseudomodel)
        os.remove(refmap_model)
    
    
        

if __name__ == '__main__':
    unittest.main()           
   
        
        
        
        
    
    
        
        
            
            
        
    
    
        
        
        
        
        
        
        
        
        
        
        
        