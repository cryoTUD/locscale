#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 17 13:28:59 2022

@author: alok
"""

import mrcfile
import os
import numpy as np
import random
from tqdm import tqdm
from locscale.include.emmer.pdb.pdb_tools import find_wilson_cutoff

folder = "/mnt/c/Users/abharadwaj1/Downloads/ForUbuntu/LocScale/tests/bfactor_correlation/refmac_iteration/"
emmap_name_list = [
    "input_map.mrc",
    "global_sharpened.mrc",
    "pseudomodel_1_refinement.mrc",
    "pseudomodel_4_refinement.mrc",
    "pseudomodel_7_refinement.mrc",
    "pseudomodel_10_refinement.mrc",
    "pseudomodel_13_refinement.mrc",
    "atomic_modeL_refined_15_iter.mrc"]

emmap_path_list = []
for emmap_name in emmap_name_list:
    emmap_path_list.append(os.path.join(folder, emmap_name))
    

correlation_map_1 = "global_sharpened.mrc"
correlation_map_2 = "pseudomodel_1_refinement.mrc"

mask_path = os.path.join(folder,"input_mask.mrc")
mask = mrcfile.open(mask_path).data

emmaps = {}
for emmap_path in emmap_path_list:
    emmap_name = os.path.basename(emmap_path)
    emmaps[emmap_name] = mrcfile.open(emmap_path).data


apix = mrcfile.open(mask_path).voxel_size.x

## Get the following from the report generated
high_frequency_cutoff = find_wilson_cutoff(mask_path=mask_path)
fsc_cutoff = 3.4
boxsize = 22#int(round(25 / apix))


def get_box(big_volume,center,size):
    return big_volume[center[2]-size//2:center[2]+size//2,center[1]-size//2:center[1]+size//2,center[0]-size//2:center[0]+size//2]

def apply_intensity_at_index(big_volume,indices,values=None):
    '''
    index : in x,y,z format
    
    '''
    for i,index in enumerate(indices):
        if values is not None:
            big_volume[index[2],index[1],index[0]] = values[i]
        else:
            big_volume[index[2],index[1],index[0]] = 1
    
    return big_volume

def get_values_at_points(big_volume,points):
    values = []
    for point in points:
        values.append(big_volume[point[2],point[1],point[0]])
    
    return values

def distance_from_center_of_box(center_of_window,shape):
    zw,yw,xw = center_of_window
    zc,yc,xc = shape[0]//2, shape[1]//2, shape[2]//2
    
    r = np.sqrt((zc-zw)**2 + (yc-yw)**2 + (xc-xw)**2)
    
    return r


#z,y,x = np.where(emmap>=0.008)
z,y,x = np.where(mask == 1)
all_points = list(zip(x,y,z))
random_centers = random.sample(all_points,15000)


bfactor_distributions = {}

analysis = {}

from locscale.include.emmer.ndimage.profile_tools import estimate_bfactor_standard, compute_radial_profile, frequency_array, plot_radial_profile
from locscale.include.emmer.ndimage.map_tools import compute_real_space_correlation

for center in tqdm(random_centers, desc="Analysing"):
    try:
        analysis[center] = []
        windows = {}
        radial_profiles = {}
        local_bfactors = {}
        
        for map_name in emmaps.keys():         
            windows[map_name] = get_box(emmaps[map_name], center, boxsize)
            radial_profiles[map_name] = compute_radial_profile(windows[map_name])
            freq = frequency_array(radial_profiles[map_name], apix)

            local_bfactors[map_name] = -1 * estimate_bfactor_standard(freq, radial_profiles[map_name],wilson_cutoff=high_frequency_cutoff, fsc_cutoff=fsc_cutoff,return_amplitude=False)
            analysis[center].append(local_bfactors[map_name])

            
    except:
        print("Error at",center)
        
import pandas as pd
import seaborn as sns
#bfactor_map = apply_intensity_at_index(np.zeros(emmap_1.shape),random_centers,bfactor_list_1)

# bin voxels with similar intensities 

num_bins = 200
import matplotlib.pyplot as plt
df = pd.DataFrame(analysis,columns=list(analysis.keys()),index=list(emmaps.keys())).T

sns.relplot(data=df,x=correlation_map_1,y=correlation_map_2)
plt.title("Bfactor correlation plots ")
plt.xlabel(" {}".format(correlation_map_1))
plt.ylabel("{}".format(correlation_map_2))
print("Correlation: {:.2f}".format(df[correlation_map_1].corr(df[correlation_map_2])))

plt.figure(2)
for map_name in emmaps.keys():
    sns.kdeplot(df[map_name])

plt.legend(list(emmaps.keys()))