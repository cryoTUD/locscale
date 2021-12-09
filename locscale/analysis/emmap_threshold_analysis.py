#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec  9 13:14:29 2021

@author: alok
"""

import mrcfile
from locscale.include.emmer.ndimage.map_quality_tools import calculate_surface_area_at_threshold, count_distinct_regions
import numpy as np
import os
from tqdm import tqdm
import matplotlib.pyplot as plt


        
    
    


folder = "/mnt/c/Users/abharadwaj1/Downloads/ForUbuntu/LocScale/tests/test_locscale/emd5778"

emmap_path = os.path.join(folder, "loc_scale_symmetrised.mrc")

emmap = mrcfile.open(emmap_path).data
apix = mrcfile.open(emmap_path).voxel_size.tolist()[0]

num_bins = 100
threshold_bins = np.linspace(0, emmap.max(), num=num_bins)
surface_area_per_threshold = {}
count_regions_per_threshold = {}
unit_surface_area_per_threshold = {}
volume_per_threshold = {}
surface_area_to_volume_threshold = {}
for threshold in tqdm(threshold_bins):
    binarised_map = (emmap>=threshold).astype(np.int_)
    sum_of_voxels = binarised_map.sum()
    surface_area_per_threshold[threshold] = calculate_surface_area_at_threshold(emmap, apix, threshold)
    count_regions_per_threshold[threshold] = count_distinct_regions(emmap, threshold)
    unit_surface_area_per_threshold[threshold] = surface_area_per_threshold[threshold] / count_regions_per_threshold[threshold]
    volume_per_threshold[threshold] = binarised_map.sum() * apix**3
    surface_area_to_volume_threshold[threshold] = surface_area_per_threshold[threshold]**3 / volume_per_threshold[threshold]**2
    
    


plt.plot(surface_area_to_volume_threshold.keys(),surface_area_to_volume_threshold.values(),'k.-')