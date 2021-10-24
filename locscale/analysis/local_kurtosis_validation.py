#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Oct 23 01:37:05 2021

@author: alok
"""

from scipy.stats import kurtosis, skew
import mrcfile
import os
import numpy as np
import random
from tqdm import tqdm
from locscale.include.emmer.ndimage.profile_tools import estimate_bfactor_standard, compute_radial_profile, frequency_array, plot_radial_profile


folder = "/mnt/c/Users/abharadwaj1/Downloads/ForUbuntu/LocScale/test_locscale/map_quality_metric"

locscale_path = os.path.join(folder,"loc_scale_oct14.mrc")

mask_path = os.path.join(folder,"pdb3j5p_mask.mrc" )

window_size = 40
wilson_cutoff = 9.6
fsc_cutoff = 3.4

def get_box(big_volume,center,size):
    return big_volume[center[2]-size//2:center[2]+size//2,center[1]-size//2:center[1]+size//2,center[0]-size//2:center[0]+size//2]

def distance_from_center_of_box(center_of_window,shape):
    zw,yw,xw = center_of_window
    zc,yc,xc = shape[0]//2, shape[1]//2, shape[2]//2
    
    r = np.sqrt((zc-zw)**2 + (yc-yw)**2 + (xc-xw)**2)
    
    return r


mask = mrcfile.open(mask_path).data
locscale_map = mrcfile.open(locscale_path).data
apix = mrcfile.open(locscale_path).voxel_size.tolist()[0]



#z,y,x = np.where(emmap>=0.008)
z,y,x = np.where(mask == 1)
all_points = list(zip(x,y,z))
random_centers = random.sample(all_points,15000)

local_analysis = {}


for center in tqdm(random_centers, desc="Validating"):
    vol = get_box(locscale_map, center, window_size)
    
    rp = compute_radial_profile(vol)
    freq = frequency_array(rp, apix)
    bfactor = estimate_bfactor_standard(freq, rp, wilson_cutoff=wilson_cutoff, fsc_cutoff=fsc_cutoff)
    sk = skew(vol.flatten())
    kurt = kurtosis(vol.flatten())
    local_analysis[center] = [bfactor, kurt, sk, vol]
    
import pandas as pd
import seaborn as sns

df = pd.DataFrame(data=local_analysis.values(), columns=['bfactor','kurtosis', 'skew','vol'])

sns.relplot(data=df, x='skew',y='kurtosis')

    


