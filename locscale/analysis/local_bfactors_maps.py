#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 21 14:43:39 2021

@author: alok
"""

import mrcfile
import os
import numpy as np
import random
from tqdm import tqdm
emdid = "0038"
folder = "/mnt/c/Users/abharadwaj1/Downloads/ForUbuntu/LocScale/tests/bfactor_correlation/emd0038/"
pseudomap_path = os.path.join(folder+"using_pseudoatomic_model_new.mrc")
atomic_map_path = os.path.join(folder+"using_atomic_model.mrc")
resmap_path = os.path.join(folder+"emd_0038_map_resmap.mrc")
mask_path = os.path.join(folder+"emd_0038_full_confidenceMap.mrc")

pseudomap = mrcfile.open(pseudomap_path).data
modmap = mrcfile.open(atomic_map_path).data
resmap = mrcfile.open(resmap_path).data
mask = mrcfile.open(mask_path).data

emmap_1 = pseudomap
emmap_2 = modmap
emmap_3 = resmap

resolution_limits = (2,6)
apix = mrcfile.open(mask_path).voxel_size.x

## Get the following from the report generated
high_frequency_cutoff = 11.58
fsc_cutoff = 3.2
boxsize = 22


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



bfactor_list_1 = []
bfactor_list_2 = []
bfactor_list_3 = []

analysis = {}

from locscale.include.emmer.ndimage.profile_tools import estimate_bfactor_standard, compute_radial_profile, frequency_array, plot_radial_profile
from locscale.include.emmer.ndimage.map_tools import compute_real_space_correlation

for center in tqdm(random_centers, desc="Analysing"):
    try:
        local_resolution = emmap_3[center[2],center[1],center[0]]
        
        vol1 = get_box(emmap_1, center, boxsize)
        vol2 = get_box(emmap_2, center, boxsize)
        vol3 = get_box(emmap_3, center, boxsize)
        
        rad_profile_1 = compute_radial_profile(vol1)
        rad_profile_2 = compute_radial_profile(vol2)
        rad_profile_3 = compute_radial_profile(vol3)
        
        freq = frequency_array(rad_profile_1, apix)
        bfactor_1, amp_1 = estimate_bfactor_standard(freq, rad_profile_1,wilson_cutoff=high_frequency_cutoff, fsc_cutoff=fsc_cutoff,return_amplitude=True)
        
        bfactor_2, amp_2 = estimate_bfactor_standard(freq, rad_profile_2,wilson_cutoff=high_frequency_cutoff, fsc_cutoff=fsc_cutoff,return_amplitude=True)
        
        bfactor_3, amp_3 = estimate_bfactor_standard(freq, rad_profile_3,wilson_cutoff=high_frequency_cutoff, fsc_cutoff=fsc_cutoff,return_amplitude=True)
        
        bfactor_list_1.append(bfactor_1)
        bfactor_list_2.append(bfactor_2)
        bfactor_list_3.append(bfactor_3)
        
        radial_distance = distance_from_center_of_box(center, emmap_3.shape)*apix
        
        rscc_profile = compute_real_space_correlation(rad_profile_1, rad_profile_2)
        
        analysis[center] = [local_resolution,bfactor_1,bfactor_2, bfactor_3, amp_1, amp_2, amp_3,radial_distance]
    except:
        print("Error at",center)
        
import pandas as pd
import seaborn as sns
#bfactor_map = apply_intensity_at_index(np.zeros(emmap_1.shape),random_centers,bfactor_list_1)

# bin voxels with similar intensities 

num_bins = 200
import matplotlib.pyplot as plt
df = pd.DataFrame(analysis,columns=list(analysis.keys()),index=['local_resolution','bfactor_1','bfactor_2','bfactor_3', 'amp_1','amp_2','amp_3','radial_distance']).T
#df['intensity_bins'] = pd.qcut(df['intensity'],q=num_bins)
#sorted_data = df.groupby(['intensity_bins']).mean()
sns.relplot(data=df,x='bfactor_1',y='bfactor_2',hue='local_resolution',hue_norm=resolution_limits)
plt.title("Correlation plots for EMD: "+emdid)
plt.xlabel("Local estimate bfactor from pseudoatomic model")
plt.ylabel("Local estimate bfactor from atomic model")
#plot_radial_profile(freq, [rad_profile_1, rad_profile_2, rad_profile_3], legends=["local pseudomap","local atomic map","local emmap"])