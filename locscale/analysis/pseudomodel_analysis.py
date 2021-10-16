#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Oct 17 00:29:54 2021

@author: alok
"""

## Pseudomodel analysis 
import os
import mrcfile

import numpy as np
import gemmi
from locscale.pseudomodel.pseudomodel_solvers import main_solver3D, main_solver_kick
from locscale.pseudomodel.pseudomodel_classes import extract_model_from_mask
from locscale.include.emmer.ndimage.map_utils import measure_mask_parameters
    
locscale_path = os.environ['LOCSCALE_PATH']

emmap_path = os.path.join(locscale_path,"tests","test_data","emd5778_map.mrc")

mask_path = os.path.join(locscale_path,"tests","test_data","emd5778_mask.mrc")

mrc = mrcfile.open(emmap_path)
emmap = mrc.data
voxelsize = float(mrc.voxel_size.x)
mask = mrcfile.open(mask_path).data
num_atoms,mask_dims = measure_mask_parameters(mask_path,verbose=True)

pseudomodel = extract_model_from_mask(mask,num_atoms,threshold=1)
    
emmap_shape = emmap.shape
unitcell = gemmi.UnitCell(emmap_shape[0]*voxelsize,emmap_shape[1]*voxelsize,emmap_shape[2]*voxelsize,90,90,90)
    

gz,gy,gx = np.gradient(emmap)
masked_grad_magnitude = mask * np.sqrt(gx**2 + gy**2 + gz**2)
max_gradient = masked_grad_magnitude.max()

g = round(100 / max_gradient)
scale_lj = 1
scale_map = 1 
friction = 10
            
save_folder = "/mnt/c/Users/abharadwaj1/Downloads/ForUbuntu/LocScale/Test/pseudomodel_analysis/"        
pseudomodel, peak_bond_length_list, map_values, cross_correlation, profiles_iterations = main_solver3D(
            emmap,gx,gy,gz,pseudomodel,g=g,friction=friction,min_dist_in_angst=1.2,voxelsize=voxelsize,dt=0.1,capmagnitude_lj=100,epsilon=1,scale_lj=scale_lj,
            capmagnitude_map=100,scale_map=scale_map,total_iterations=50, compute_map=True,emmap_path=None,mask_path=None,returnPointsOnly=False,
            integration='verlet',verbose=True, save_path=save_folder)


from locscale.include.emmer.ndimage.profile_tools import frequency_array, plot_radial_profile, compute_radial_profile
import matplotlib.pyplot as plt
rp_emmap = compute_radial_profile(emmap)
freq = frequency_array(rp_emmap, apix=float(voxelsize))
plt.figure(1)
plot_radial_profile(freq, profiles_iterations)
plt.savefig("profile_iterations_correlation.jpg")

plt.figure(2)
plt.plot(cross_correlation)
plt.xlabel("Pseudomodel iterations")
plt.ylabel("RSCC with EM map")

plt.figure(3)
map_val, sd = zip(*map_values)
plt.errorbar(np.linspace(1,51,50), map_val, sd)
plt.xlabel("Pseudomodel iterations")
plt.ylabel("Average map value")



