#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar  4 15:18:27 2022

@author: alok
"""

import mrcfile
import seaborn as sns
import numpy as np
import os
import matplotlib.pyplot as plt

def get_local_resolution_data(local_resolution_map, mask, mask_threshold):
    from locscale.include.emmer.ndimage.map_utils import parse_input
    
    local_resolution_map = parse_input(local_resolution_map)
    mask = parse_input(mask)
    
    binarised_LR_map = (mask>=mask_threshold).astype(np.int_)
    flattend_array = (binarised_LR_map * local_resolution_map).flatten()
    nonzero_array = flattend_array[flattend_array>0]
    
    return nonzero_array
    
def get_values_at_points(big_volume,points):
    values = []
    for point in points:
        values.append(big_volume[point[2],point[1],point[0]])
    
    return values

folder = "/mnt/c/Users/abharadwaj1/Downloads/ForUbuntu/Resolution_testing/local_resolution_analysis"

unsharpened_map_path = os.path.join(folder, "emd3061_unsharpened_halfmaps_fullmap.mrc")
mask_path = os.path.join(folder, "emd_3061_confidence_final.mrc")
deepEmhancer_localRes_path = os.path.join(folder, "deepemhancer_halfmap_1_localResolutions.mrc")
unsharpened_localRes_path = os.path.join(folder, "emd3061_unsharpened_localresolutions.mrc")
locscale_localRes_path = os.path.join(folder, "emd3061_locscale_sharpened_localresolutions.mrc")
nn_localRes_path = os.path.join(folder, "emd3061_nn_pred_localresolution.mrc")

no_lipid_threshold = 0.021

## Start calculation

emmap = mrcfile.open(unsharpened_map_path).data
binarised_mask_no_lipid = (emmap>no_lipid_threshold).astype(np.int_)

unsharp_locres_array = get_local_resolution_data(unsharpened_localRes_path, binarised_mask_no_lipid, mask_threshold=1)
locscale_sharp_locres_array = get_local_resolution_data(locscale_localRes_path, binarised_mask_no_lipid, mask_threshold=1)
nn_pred_locres_array = get_local_resolution_data(nn_localRes_path, binarised_mask_no_lipid, mask_threshold=1)
deepEmhancer_locres_array = get_local_resolution_data(deepEmhancer_localRes_path, binarised_mask_no_lipid, mask_threshold=1)

sns.kdeplot(unsharp_locres_array), sns.kdeplot(locscale_sharp_locres_array), sns.kdeplot(nn_pred_locres_array), sns.kdeplot(deepEmhancer_locres_array)
plt.legend(["Unsharpened halfmaps", "Locscale sharpened halfmaps", "Emmernet halfmaps","DeepEmhancer halfmaps"])
plt.xlabel("Local resolution ($\AA$)")


#%% Finding correlations in the data
import random
from locscale.include.emmer.ndimage.map_utils import get_all_voxels_inside_mask
import pandas as pd

all_inside_mask = get_all_voxels_inside_mask(binarised_mask_no_lipid, 1)

sample_size = 10000

random_centers = random.sample(all_inside_mask, sample_size)

local_resolutions_unsharp = get_values_at_points(mrcfile.open(unsharpened_localRes_path).data, random_centers)
local_resolutions_locscale = get_values_at_points(mrcfile.open(locscale_localRes_path).data, random_centers)
local_resolutions_nn = get_values_at_points(mrcfile.open(nn_localRes_path).data, random_centers)

df = pd.DataFrame(data=[local_resolutions_unsharp, local_resolutions_locscale, local_resolutions_nn], index=["Unsharp","Locscale","NN"]).T






