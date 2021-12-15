#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec 10 19:48:42 2021

@author: alok
"""

from scipy.signal import argrelextrema
import mrcfile
import numpy as np
import csv
import pickle
import os
import matplotlib.pyplot as plt

folder = "/mnt/c/Users/abharadwaj1/Downloads/ForUbuntu/LocScale/Test/threshold_analysis"
pickle_filename = os.path.join(folder, "threshold_analysis_combined.pickle")

def plot_dictionary(dictionary):
    x = np.array(list(dictionary.keys()))
    y = np.array(list(dictionary.values()))
    plt.plot(x, y, 'k.-')
    
with open(pickle_filename, "rb+") as file:
    threshold_analysis = pickle.load(file)

inv_compactness_curves = {}    
inv_compactness_curves_native = {}
local_minima = {}
emdbs_with_no_minima = []
for emdb_pdb in threshold_analysis.keys():
    unsharp_result = threshold_analysis[emdb_pdb]['unsharp_analysis']
    thresholds_native = np.array(list(unsharp_result.keys()))
    inv_compactness_native = np.array(list(unsharp_result.values()))
    local_minima[emdb_pdb] = argrelextrema(inv_compactness_native, np.less)
    first_minima = local_minima[emdb_pdb][0]
    if len(first_minima) > 0:
        value_at_first_minima = inv_compactness_native[first_minima[0]]
        if value_at_first_minima > 25:
            ignore_first = first_minima[0]
        else:
            ignore_first = 0
            emdbs_with_no_minima.append(emdb_pdb)
    else:
        ignore_first = 0
        emdbs_with_no_minima.append(emdb_pdb)
        
    
    max_after_initial_drop = inv_compactness_native[ignore_first:].max()
    
    inv_compactness_norm = inv_compactness_native/max_after_initial_drop
    
    inv_compactness_curves[emdb_pdb] = inv_compactness_norm[ignore_first:]
    inv_compactness_curves_native[emdb_pdb] = {
        'x':thresholds_native,
        'y':inv_compactness_native}
    

#%%
import scipy.cluster.hierarchy as hac
from scipy.cluster.hierarchy import fcluster
import pandas as pd
from scipy.spatial.distance import pdist

profiles_df = pd.DataFrame(inv_compactness_curves.values(), index=inv_compactness_curves.keys()).T

distance_matrix = profiles_df.corr()

condensed_distance_matrix = pdist(distance_matrix)
Z = hac.linkage(condensed_distance_matrix, method='complete',metric='euclidean')

cluster = fcluster(Z, 4, criterion='maxclust')
    
#%%

cluster_curves= {}
cluster_keys = {}
max_cluster_groups = cluster.max()

## initialise cluster curves
for label in range(1, max_cluster_groups+1):
    cluster_curves[label] = []
    cluster_keys[label] = []
## Append curves    
for i, current_emdb in enumerate(inv_compactness_curves.keys()):
    current_label = cluster[i]
    current_curve = inv_compactness_curves[current_emdb]
    for label in range(1, max_cluster_groups+1):
        if current_label == label:
            cluster_curves[label].append(current_curve)
            cluster_keys[label].append(current_emdb)

#%% 
## Plot curves

for label in range(1, max_cluster_groups+1):
    fig = plt.figure()
    plt.title("Threshold curves for label: {}".format(label))
    for plot_curve in cluster_curves[label]:
        plt.plot(plot_curve, 'k', alpha=0.2)
    
    

    
            


    
    