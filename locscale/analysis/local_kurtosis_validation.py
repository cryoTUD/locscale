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
from locscale.include.emmer.pdb.pdb_tools import find_wilson_cutoff

folder = "/mnt/c/Users/abharadwaj1/Downloads/ForUbuntu/LocScale/tests/map_quality/xray_maps"

locscale_path = os.path.join(folder,"emd_21109_map_centered.mrc")
emmap_path = os.path.join(folder, "emd_21109_map_centered.mrc")
mask_path = os.path.join(folder,"emd_21109_mask.mrc" )

window_size = 40
wilson_cutoff = find_wilson_cutoff(mask_path=mask_path)
fsc_cutoff = 1.5

def get_box(big_volume,center,size):
    return big_volume[center[2]-size//2:center[2]+size//2,center[1]-size//2:center[1]+size//2,center[0]-size//2:center[0]+size//2]

def distance_from_center_of_box(center_of_window,shape):
    zw,yw,xw = center_of_window
    zc,yc,xc = shape[0]//2, shape[1]//2, shape[2]//2
    
    r = np.sqrt((zc-zw)**2 + (yc-yw)**2 + (xc-xw)**2)
    
    return r


mask = mrcfile.open(mask_path).data
locscale_map = mrcfile.open(locscale_path).data
emmap = mrcfile.open(emmap_path).data

apix = mrcfile.open(locscale_path).voxel_size.tolist()[0]



#z,y,x = np.where(emmap>=0.008)
z,y,x = np.where(mask == 1)
all_points = list(zip(x,y,z))
random_centers = random.sample(all_points,15000)

local_analysis = {}

volumes = {}
for center in tqdm(random_centers, desc="Validating"):
    try:
        distance_to_center = distance_from_center_of_box(center, locscale_map.shape)
       
        window_locscale = get_box(locscale_map, center, window_size)
        window_emmap = get_box(emmap, center, window_size)
        
        ## calculate rp
        rp_locscale = compute_radial_profile(window_locscale)
        rp_emmap = compute_radial_profile(window_emmap)
        freq = frequency_array(rp_locscale, apix)
        bfactor_locscale = estimate_bfactor_standard(freq, rp_locscale, wilson_cutoff=wilson_cutoff, fsc_cutoff=fsc_cutoff)
        bfactor_emmap = estimate_bfactor_standard(freq, rp_emmap, wilson_cutoff=wilson_cutoff, fsc_cutoff=fsc_cutoff)
        
        ## histogram metrics
        
        skew_locsale = skew(window_locscale.flatten())
        kurtosis_locscale = kurtosis(window_locscale.flatten())
        skew_emmap = skew(window_emmap.flatten())
        kurtosis_emmap = kurtosis(window_emmap.flatten())
        
        
        
        local_analysis[center] = [bfactor_locscale, bfactor_emmap, kurtosis_locscale, skew_locsale, kurtosis_emmap, skew_emmap, tuple(center), distance_to_center]
        volumes[center] = window_locscale
    except:
        continue


 
#%%
import pandas as pd
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
from datetime import datetime
import mplcursors
import seaborn as sns

df = pd.DataFrame(data=local_analysis.values(), columns=['bfactor_locscale', 'bfactor_emmap','kurtosis_locscale', 'skew_locscale',
                                                         'kurtosis_emmap', 'skew_emmap','center', 'radius'])

def binomial_quadratic(x):
    return 1+x**2

def general_quadratic(x,a,b,c):
    return a + b*x + c* x**2

def r2(y_fit, y_data):
    y_mean = y_data.mean()
    residual_squares = (y_data-y_fit)**2
    variance = (y_data-y_mean)**2
    
    residual_sum_of_squares = residual_squares.sum()
    sum_of_variance = variance.sum()
    
    r_squared = 1 - residual_sum_of_squares/sum_of_variance
    
    return r_squared

df.sort_values(by=['skew_locscale'])

def clickable_plot(data, x_col, y_col, click_val):
    df = data.copy()

    labels = list(df[click_val])
    y_fit = binomial_quadratic(df[x_col].sort_values())
    r_squared = r2(y_fit, df[y_col])
    
    fig, ax = plt.subplots(1, figsize=(8,6))
    scatter_plot=ax.scatter(df[x_col], df[y_col])
    cursor = mplcursors.cursor(scatter_plot, hover=False)
    plt.xlabel(x_col)
    plt.ylabel(y_col)
    
    @cursor.connect('add')
    def on_add(sel):
        x_col_val = sel.target[0]
        y_col_val = sel.target[1]
        click_value = df.iloc[sel.target.index][click_val]
        print("Center: {} | {}: {:.2f} & {}: {:.2f}  ".format(click_value, x_col, x_col_val, y_col, y_col_val))
        sel.annotation.set(text=click_value)
        
    
    plt.show()


#%%

'''
def click(event):
    if event.inaxes == ax:
       
        cont, ind = scatter_plot.contains(event)
        
        if cont:
            annot.xy = (event.xdata, event.ydata)
            print(event.xdata, event.ydata)
            annot.set_text("hello")
            annot.set_visible(True)
        else:
            annot.set_visible(False)
            
            
            
'''    
