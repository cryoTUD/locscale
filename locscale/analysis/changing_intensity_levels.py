#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 15 10:37:24 2021

@author: alok
"""

import mrcfile
import gemmi
import numpy as np
from emmer.ndimage.profile_tools import compute_radial_profile, plot_radial_profile, frequency_array
from emmer.ndimage.map_utils import save_as_mrc
from locscale.pseudomodel.pseudomodel_headers import normalise_intensity_levels

folder = "/mnt/c/Users/abharadwaj1/Downloads/ForUbuntu/LocScale/locscale_check/pdb2map_cursed/changing_intensity_levels/"
test_emmap_path = folder+"emd5778_unfiltered.mrc"
correct_emmap_path = folder+"emd5778_tutorial.mrc"
gemmi_refmap_path = folder + "pdb3j5p_using_pdb2map.mrc"

test_emmap = mrcfile.open(test_emmap_path).data
correct_emmap = mrcfile.open(correct_emmap_path).data
gemmi_refmap = mrcfile.open(gemmi_refmap_path).data


## Generate datasets: 
    #1 Normalise test_emmap_path to gemmi_refmap levels
    #2 Normalise correct_emmap_path to gemmi_refmap levels
    #2a

    #3 Iteratively normalise gemmi_refmap to test_emmap level in ten steps

gemmi_levels = np.array([gemmi_refmap.min(),gemmi_refmap.max()])
correct_emmap_levels = np.array([correct_emmap.min(),correct_emmap.max()])
test_emmap_levels = np.array([test_emmap.min(),test_emmap.max()])

normalised_test_emmap = normalise_intensity_levels(
    test_emmap, gemmi_levels)

normalised_correct_emmap = normalise_intensity_levels(
    correct_emmap, to_levels=gemmi_levels)

normalised_correct_emmap_high = normalise_intensity_levels(
    correct_emmap, to_levels=correct_emmap_levels*1e6)


rp_emmap = compute_radial_profile(test_emmap)
rp_correct_emmap = compute_radial_profile(correct_emmap)
rp_gemmi_map = compute_radial_profile(gemmi_refmap)
rp_normalised_test_emmap = compute_radial_profile(normalised_test_emmap)
rp_normalised_correct_emmap = compute_radial_profile(normalised_correct_emmap)
rp_normalised_correct_emmap_high = compute_radial_profile(normalised_correct_emmap_high)


scale_factor = (rp_emmap / rp_gemmi_map).mean()
normalised_gemmi_maps = {}
rp_normalised_gemmi_maps = {}

for i in np.linspace(0,1,10):    
    normalised_gemmi_map = normalise_intensity_levels(
        gemmi_refmap, to_levels=gemmi_levels*scale_factor*i)
    normalised_gemmi_maps[i] = normalised_gemmi_map

    rp_normalised_gemmi_maps[i] = compute_radial_profile(normalised_gemmi_map)

freq = frequency_array(rp_emmap, 1.2156)

plot_radial_profile(freq, list(rp_normalised_gemmi_maps.values())[:5])


save_as_mrc(normalised_test_emmap, apix=1.2156, output_filename=folder+"normalised_test_emmap.mrc")
save_as_mrc(normalised_correct_emmap, apix=1.2156, output_filename=folder+"normalised_correct_emmap.mrc")
save_as_mrc(normalised_correct_emmap_high, apix=1.2156, output_filename=folder+"normalised_correct_emmap_high.mrc")

for i in np.linspace(0,1,10):
    filename = folder+"normalised_gemmi_maps_"+str(round(i*100))+".mrc"
    
    
    save_as_mrc(normalised_gemmi_maps[i], apix=1.2156, output_filename=filename)


    




    
    
                     


