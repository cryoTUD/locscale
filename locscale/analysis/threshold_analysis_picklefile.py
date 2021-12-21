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

folder = "/mnt/c/Users/abharadwaj1/Downloads/ForUbuntu/LocScale/Test/threshold_analysis/"
pickle_filename = os.path.join(folder, "threshold_analysis_combined_DC_metric.pickle")

def normalise(x):
    return x/x.max()
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
    thresholds_native = np.array(list(unsharp_result['threshold']))
    avg_length = np.array(list(unsharp_result['avg_lengths_array']))
    total_surface_area = unsharp_result['total_surface_area_array']
    total_volume= unsharp_result['total_volume_array']
    detail = total_surface_area/total_volume
    connection = avg_length
    
    dc_metric_not_normal = detail * connection
    
    local_minima[emdb_pdb] = argrelextrema(dc_metric_not_normal, np.less)
    first_minima = local_minima[emdb_pdb][0]
    '''
    if len(first_minima) > 0:
        value_at_first_minima = dc_metric_not_normal[first_minima[0]]
        if value_at_first_minima > 25:
            ignore_first = first_minima[0]
        else:
            ignore_first = 0
            emdbs_with_no_minima.append(emdb_pdb)
    else:
        ignore_first = 0
        emdbs_with_no_minima.append(emdb_pdb)
     '''
     
    
    #max_after_initial_drop = inv_compactness_native[ignore_first:].max()
    
    #inv_compactness_norm = inv_compactness_native/max_after_initial_drop
    ignore_first = first_minima[0]
    from locscale.include.emmer.ndimage.profile_tools import resample_1d
    dc_metric = dc_metric_not_normal
    dc_metric[:ignore_first] = 0
    x = np.arange(ignore_first, 100, 1)
    y = dc_metric[ignore_first:]
    
    dc_metric_resampled = resample_1d(x, y, 100, )[1]
    dc_metric_resampled_norm = normalise(dc_metric_resampled)
    inv_compactness_curves[emdb_pdb] = dc_metric_resampled_norm
    inv_compactness_curves_native[emdb_pdb] = {
        'x':thresholds_native,
        'y':dc_metric_not_normal}
    

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
    
    
#%%
def find_optimal_threshold(curves_dict, emdb_pdb):
    x = curves_dict[emdb_pdb]['x']
    y = curves_dict[emdb_pdb]['y']
    ignore_first = 5#local_minima[emdb_pdb][0][0]
    
    y_max = y[ignore_first:].max()
    index = np.where(y==y_max)
    
    fig, ax = plt.subplots()
    ax.plot(x[ignore_first:], y[ignore_first:], 'k.-')
    threshold = x[index]
    ax.set_title("Threshold for {} is {}".format(emdb_pdb, threshold))
    return fig

def mean_top_pc(x, threshold=0.01):
    n = len(x)
    one_pc_n = int(n*threshold)
    if one_pc_n < 1:
        one_pc_n = 1
    top_pc_array = x[x.argsort()[-1*one_pc_n:][::-1]]
    
    mean_top_pc = top_pc_array.mean()
    
    return mean_top_pc
from statistics import median
def plot_threshold_curves(threshold_analysis, emdb_pdb, curve_type="unsharp_analysis", noise_threshold=0):
    import matplotlib.pyplot as plt
    x = threshold_analysis[emdb_pdb][curve_type]['threshold']
    lengths_array = threshold_analysis[emdb_pdb][curve_type]['lengths_array']
    max_length = np.array([mean_top_pc(x,0.5) for x in lengths_array])
    tsa = threshold_analysis[emdb_pdb][curve_type]['total_surface_area_array']
    tv = threshold_analysis[emdb_pdb][curve_type]['total_volume_array']
    num_regions = threshold_analysis[emdb_pdb][curve_type]['num_regions_array']
    x = x[noise_threshold:]
    tsa = tsa[noise_threshold:]
    tv = tv[noise_threshold:]
    max_length = max_length[noise_threshold:]
    num_regions = num_regions[noise_threshold:]
    
    detail = tsa/tv
    connection = normalise(max_length)
    dc_metric = detail * connection
    
    fig, ax = plt.subplots(nrows=3)
    fig.suptitle("Threshold analysis: {}, {}".format(emdb_pdb, curve_type))
    ax[0].plot(x, detail, 'b.-')
    ax[0].set_ylabel('Detail (b)')
    ax[0].get_xaxis().set_visible(False)
    #ax[0].set_ylim([0, 1])
    
    ax2 = ax[0].twinx()
    ax2.plot(x, connection, 'r.-')
    ax2.set_ylabel("Connection (r)")
    ax2.get_xaxis().set_visible(False)
   # ax2.set_ylim([0, 1])
    
    ax[1].plot(x, np.log(num_regions), 'k.-')
    ax[1].set_ylabel('log(#regions)')
    ax[1].get_xaxis().set_visible(False)
    
    ax[2].plot(x, dc_metric, 'k.-')
    ax[2].set_ylabel("DC metric")
    #ax[2].set_ylim([0, 1])
    
    plt.show()
    plt.tight_layout()
    
    return fig
    

def save_curves_with_threshold(curves_dict, filename="curves.pdf"):
    from matplotlib.backends.backend_pdf import PdfPages
    pfile = PdfPages(filename)
    import random
    for emdb in curves_dict.keys():
        
        noise = local_minima[emdb][0][0]
        
        
        fig_unsharp = plot_threshold_curves(curves_dict, emdb, noise_threshold=noise)
        #fig_sharp = plot_threshold_curves(curves_dict, emdb, curve_type="sharp_analysis", noise_threshold=noise)
        pfile.savefig(fig_unsharp)
        #pfile.savefig(fig_sharp)
    
    pfile.close()


    
    