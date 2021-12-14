#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 25 16:24:30 2021

@author: alok
"""

import mrcfile
import numpy as np
import csv
import pickle
import os

def mesh_surface_area(data, threshold, apix):
    from skimage import measure
   
    verts, faces,_,_ = measure.marching_cubes(data, threshold)
    surface_area = measure.mesh_surface_area(verts, faces) * apix**2
    return surface_area

def find_volume_matching_threshold(emmap, reference_volume, apix, num_bins=100):
    threshold_bins = np.linspace(0, emmap.max(), num=num_bins)
    
    for threshold in threshold_bins:
        binarised_map = (emmap>=threshold).astype(np.int_)
        sum_of_voxels = binarised_map.sum()
        volume_real_units = sum_of_voxels * apix**3
        if volume_real_units <= reference_volume:
            matching_threshold = threshold
            break
    
    
    return matching_threshold

def calculate_blob_surface_statistics(emmap_path, mask_path, mask_emmap=True):
    from skimage.morphology import skeletonize
    def count_distinct_regions_inner(emmap, reference_threshold):
        from skimage import measure
        
        binarised_emmap = (emmap>reference_threshold).astype(np.int_)
        labels, num_regions = measure.label(binarised_emmap, background=0, return_num=True)
        
        return labels, num_regions
    
    import mrcfile
    print("Calculating blob surface statistics for: {} using mask {}".format(emmap_path.split("/")[-1], mask_path.split("/")[-1]))
    emmap = mrcfile.open(emmap_path).data
    apix = tuple(mrcfile.open(emmap_path).voxel_size.tolist())[0]
    apix = apix
    origin = mrcfile.open(emmap_path).header.origin.tolist()
    
    mask = mrcfile.open(mask_path).data
    mask = (mask==1).astype(np.int_)
    
    if mask_emmap:
        emmap = emmap * mask
    
    mask_volume = mask.sum() * apix**3
    reference_mask_volume = mask_volume * 0.2  ## Thresholded at 20% of molecular volume  
    
    print("Finding reference threshold corresponding to 20% of molecular volume determined from mask = {} ang cubed ".format(reference_mask_volume))
    reference_threshold = find_volume_matching_threshold(emmap, reference_mask_volume, apix)
    print("Reference threshold found to be {:.2f}".format(reference_threshold))
    
    labels, num_regions = count_distinct_regions_inner(emmap, reference_threshold)
    
    surface_area_per_region = {}
    volume_per_region = {}
    surface_area_to_volume = {}
    length_per_region = {}
    
    print("Number of regions: ", num_regions)
    
    for region in np.unique(labels):
        binarised_label = (labels==region).astype(np.int_)
        surface_area_per_region[region] = mesh_surface_area(binarised_label, 0.999999, apix)
        volume_per_region[region] = binarised_label.sum() * apix**3
        surface_area_to_volume[region] = surface_area_per_region[region] /volume_per_region[region]
        binarised_for_skeletonize = binarised_label
        binarised_for_skeletonize[binarised_for_skeletonize==0] = -1
        skeletonised_label = skeletonize(binarised_for_skeletonize)
        length_per_region[region] = skeletonised_label.sum() * apix
        
    
    binarised_emmap = (emmap>reference_threshold).astype(np.int_)
    total_surface_area = mesh_surface_area(binarised_emmap, 0.999999, apix)
    
    surface_area_to_volume_array = np.array(list(surface_area_to_volume.values()))
    surface_area_array = np.array(list(surface_area_per_region.values()))
    volume_array = np.array(list(volume_per_region))
    
    blob_statistics = {
        'total_surface_area':total_surface_area,
        'num_regions':num_regions,
        'total_volume':mask_volume,
        'surface_area_to_volume_array':surface_area_to_volume_array,
        'surface_area_array':surface_area_array,
        'volume_array':volume_array,
        'length_region':length_per_region}
    
    return blob_statistics

def calculate_surface_statistics_threshold(emmap_path):
    from tqdm import tqdm
    from skimage.morphology import skeletonize
    def count_distinct_regions_inner(emmap, reference_threshold):
        from skimage import measure
        
        binarised_emmap = (emmap>reference_threshold).astype(np.int_)
        labels, num_regions = measure.label(binarised_emmap, background=0, return_num=True)
        
        return labels, num_regions
    
    import mrcfile
    print("Calculating surface statistics for: {} at different threshold".format(emmap_path.split("/")[-1]))
    emmap = mrcfile.open(emmap_path).data
    apix = tuple(mrcfile.open(emmap_path).voxel_size.tolist())[0]
    apix = apix
    origin = mrcfile.open(emmap_path).header.origin.tolist()
    
    
    num_bins = 100
    threshold_bins = np.linspace(0, emmap.max()*0.99, num=num_bins)
    surface_stats_threshold = {}
    for reference_threshold in tqdm(threshold_bins):
        labels, num_regions = count_distinct_regions_inner(emmap, reference_threshold)
        
                   
        
        binarised_emmap = (emmap>reference_threshold).astype(np.int_)
        total_surface_area = mesh_surface_area(binarised_emmap, 0.999999, apix)
        total_volume = binarised_emmap.sum() * apix**3
        
        unit_surface_area = total_surface_area / num_regions
        inverse_compactness = total_surface_area**3 / total_volume**2
       
        
        surface_stat = {
            'total_surface_area':total_surface_area,
            'num_regions':num_regions,
            'total_volume':total_volume,
            'unit_surface_area':unit_surface_area,
            'inverse_compactness':inverse_compactness}
        
        surface_stats_threshold[reference_threshold] = surface_stat
    
    total_surface_area_array = np.array([x['total_surface_area'] for x in surface_stats_threshold.values()]) 
    num_regions_array = np.array([x['num_regions'] for x in surface_stats_threshold.values()]) 
    total_volume_array = np.array([x['total_volume'] for x in surface_stats_threshold.values()]) 
    unit_surface_area_array = np.array([x['unit_surface_area'] for x in surface_stats_threshold.values()]) 
    inverse_compactness_array = np.array([x['inverse_compactness'] for x in surface_stats_threshold.values()]) 
    threshold = np.array([x for x in surface_stats_threshold.keys()]) 
    
    surface_statistics_compiled = {
        'total_surface_area_array':total_surface_area_array,
        'num_regions_array':num_regions_array,
        'total_volume_array':total_volume_array,
        'unit_surface_area_array':unit_surface_area_array,
        'inverse_compactness_array':inverse_compactness_array,
        'threshold':threshold}
        
    return surface_statistics_compiled
def threshold_analysis(emmap_path):
    from tqdm import tqdm
    from locscale.include.emmer.ndimage.map_quality_tools import calculate_surface_area_at_threshold, count_distinct_regions

    emmap = mrcfile.open(emmap_path).data
    apix = mrcfile.open(emmap_path).voxel_size.tolist()[0]
    
    num_bins = 100
    threshold_bins = np.linspace(0, emmap.max(), num=num_bins)
    
    surface_area_to_volume_threshold = {}
    for threshold in threshold_bins:
        binarised_map = (emmap>=threshold).astype(np.int_)
        sum_of_voxels = binarised_map.sum()
        surface_area_per_threshold = calculate_surface_area_at_threshold(emmap, apix, threshold)
        
        volume_per_threshold = binarised_map.sum() * apix**3
        surface_area_to_volume_threshold[threshold] = surface_area_per_threshold**3 / volume_per_threshold**2
        
    return surface_area_to_volume_threshold

def plot_threshold_analysis_data(x,y, emdb_filename, crop_from=4):
    from matplotlib.offsetbox import AnchoredText
    import matplotlib.pyplot as plt
    
    max_y_value = y[crop_from:].max()  ## Excluding first few threshold due to noise
    print(np.where(y==max_y_value)[0][0])
    threshold = x[np.where(y==max_y_value)[0][0]]
            
    fig, ax1 = plt.subplots(1,1)
    ax1.plot(x[crop_from:], y[crop_from:],'k.-')
    ax1.set_xlabel("Threshold")
    ax1.set_ylabel("SA$^3$ / V$^2$")
    ax1.set_title("Surface area to volume threshold for {}".format(emdb_filename))
    text = "Optimal threshold at: {}".format(round(threshold, 4))
    anchored_text=AnchoredText(text, loc=3)
    ax1.add_artist(anchored_text)
    
    return threshold, fig

def plot_surface(emmap, threshold):
    
    from matplotlib.offsetbox import AnchoredText
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d.art3d import Poly3DCollection
    from skimage import measure  
    from scipy.ndimage import rotate
    
    rotated_emmap = rotate(emmap, 90, axes=(0,2))
    try:
        binarised_emmap = (rotated_emmap>=threshold).astype(np.int_)
        x,y,z = np.where(binarised_emmap==1)
        verts, faces, _, __ = measure.marching_cubes(binarised_emmap, 0.99999)
        mesh = Poly3DCollection(verts[faces])
        mesh.set_alpha(0.1)
        mesh.set_facecolor('0.2')
        fig = plt.figure()
        
        ax = fig.add_subplot(111, projection='3d')
        ax.set_axis_off()
        ax.add_collection3d(mesh)
        text = "Optimal threshold at: {}".format(round(threshold, 4))
        anchored_text=AnchoredText(text, loc=2)
        ax.add_artist(anchored_text)
        ax.set_xlim(x.min(), x.max())  
        ax.set_ylim(y.min(), y.max())  
        ax.set_zlim(z.min(), z.max())  
        

        plt.tight_layout()
        return fig
    except Exception as e:
        
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        text = "Oops! something went wrong"
        anchored_text=AnchoredText(text, loc=2)
        ax.add_artist(anchored_text)
        print(e)
        return fig
    


def get_threshold_analysis_data(parent_folder, EMDB_PDB_ids,output_filename, process_id=None):
    #from matplotlib.backends.backend_pdf import PdfPages
    import matplotlib.pyplot as plt
    from matplotlib.offsetbox import AnchoredText
    
    #pdf_file = PdfPages(output_filename)
    
        
    unsharpened_map_path_prefix = "EMDBmaps/"
    sharpened_map_path_prefix = "EMDBmaps/sharpened_map/"
    pdb_path_prefix = "PDBgemmis/"
    
    threshold_statistics = {}
    
    pickle_output_file = os.path.join(parent_folder, "{}.pickle".format(output_filename))
    #csv_output_file = os.path.join(parent_folder, "{}_csv.csv".format(output_filename))
    
    for emdb_pdb in EMDB_PDB_ids:
        print("#################################################################################################### \n")
        print("#################################################################################################### \n")
        #csv_writer = open(csv_output_file, mode="a")
        print(emdb_pdb)
        try:
                
            emdb_id = emdb_pdb.split("_")[0]
            pdb_id = emdb_pdb.split("_")[1]
    
            
            emdb_filename = "emd_{}_unsharpened.map".format(emdb_id)
            mask_filename = "emd_{}_confidence.map".format(emdb_id)
            pdb_filename = "PDB_{}.pdb".format(pdb_id)
            MD_locscale_filename = "emd_{}_MD_sharpened_refit_10_blur_20.map".format(emdb_id)
            
            unsharpened_map_path = os.path.join(parent_folder, emdb_pdb, unsharpened_map_path_prefix, emdb_filename)
            mask_path = os.path.join(parent_folder, emdb_pdb, unsharpened_map_path_prefix, mask_filename)
            md_locscale_path = os.path.join(parent_folder, emdb_pdb, sharpened_map_path_prefix, MD_locscale_filename)
            
            emmap_unsharpened = mrcfile.open(unsharpened_map_path).data
            emmap_sharpened = mrcfile.open(md_locscale_path).data
                    
            surface_area_to_volume_ratio_unsharp = threshold_analysis(unsharpened_map_path)
            surface_area_to_volume_ratio_sharp = threshold_analysis(md_locscale_path)
            
            x_unsharp  = np.array(list(surface_area_to_volume_ratio_unsharp.keys()))
            y_unsharp = np.array(list(surface_area_to_volume_ratio_unsharp.values()))
            
      
            x_sharp = np.array(list(surface_area_to_volume_ratio_sharp.keys()))
            y_sharp = np.array(list(surface_area_to_volume_ratio_sharp.values()))
            
            
            unsharp_threshold, unsharp_plot_xy = plot_threshold_analysis_data(x_unsharp, y_unsharp, 
                                                                              emdb_filename)
            
            sharp_threshold, sharp_plot_xy = plot_threshold_analysis_data(x_sharp, y_sharp, 
                                                                              MD_locscale_filename)
            
            #unsharp_surface_fig = plot_surface(emmap_unsharpened, unsharp_threshold)
            #sharp_surface_fig = plot_surface(emmap_sharpened, sharp_threshold)
            
            unsharp_plot_xy.savefig("{}_1_unsharpened_threshold_analysis.png".format(emdb_pdb))
            sharp_plot_xy.savefig("{}_2_sharpened_threshold_analysis.png".format(emdb_pdb))
            #unsharp_surface_fig.savefig("{}_3_unsharpened_surface_plot.png".format(emdb_pdb))
            #sharp_surface_fig.savefig("{}_4_sharpened_surface_plot.png".format(emdb_pdb))
            print("saved plots for {}".format(emdb_pdb))
            threshold_statistics[emdb_pdb] = {'unsharp_analysis':surface_area_to_volume_ratio_unsharp, 'sharp_analysis':surface_area_to_volume_ratio_sharp}


        except Exception as e:
           # fig_error = plt.figure()
           # ax = fig_error.add_subplot(111, projection='3d')
           # text = "Oops! something went wrong {}".format(emdb_pdb)
           # anchored_text=AnchoredText(text, loc=2)
           # ax.add_artist(anchored_text)
            print("Error for {}".format(emdb_pdb))
            print(e)
    with open(pickle_output_file, "wb") as picklefile:
        pickle.dump(threshold_statistics, picklefile)
           
            
            
    
def get_data(parent_folder, EMDB_PDB_ids, fsc_resolutions_list,output_filename, process_id=None):
    
    fsc_resolution = {}
    for data in fsc_resolutions_list:
        pdb_id = data.split("_")[0]
        resolution = float(data.split("_")[1])
        fsc_resolution[pdb_id] = resolution
        
    unsharpened_map_path_prefix = "EMDBmaps/"
    sharpened_map_path_prefix = "EMDBmaps/sharpened_map/"
    pdb_path_prefix = "PDBgemmis/"
    
    map_blob_statistics = {}
    
    pickle_output_file = os.path.join(parent_folder, "{}.pickle".format(output_filename))
    csv_output_file = os.path.join(parent_folder, "{}_csv.csv".format(output_filename))
    
    for emdb_pdb in EMDB_PDB_ids:
        print("#################################################################################################### \n")
        print("#################################################################################################### \n")
        csv_writer = open(csv_output_file, mode="a")
        print(emdb_pdb)
        try:
                
            emdb_id = emdb_pdb.split("_")[0]
            pdb_id = emdb_pdb.split("_")[1]
    
            
            emdb_filename = "emd_{}_unsharpened.map".format(emdb_id)
            mask_filename = "emd_{}_confidence.map".format(emdb_id)
            pdb_filename = "PDB_{}.pdb".format(pdb_id)
            MD_locscale_filename = "emd_{}_MD_sharpened_refit_10_blur_20.map".format(emdb_id)
            
            unsharpened_map_path = os.path.join(parent_folder, emdb_pdb, unsharpened_map_path_prefix, emdb_filename)
            mask_path = os.path.join(parent_folder, emdb_pdb, unsharpened_map_path_prefix, mask_filename)
            md_locscale_path = os.path.join(parent_folder, emdb_pdb, sharpened_map_path_prefix, MD_locscale_filename)
                    
    
            unsharpened_map_blob_statistics = calculate_blob_surface_statistics(unsharpened_map_path, mask_path)
            sharpened_blob_statistics = calculate_blob_surface_statistics(md_locscale_path, mask_path)
            
            
            map_blob_statistics[emdb_pdb] = {
                'unsharpened_map_blob_statistics':unsharpened_map_blob_statistics, 
                'sharpened_blob_statistics':sharpened_blob_statistics,
                }
            
            for key in map_blob_statistics[emdb_pdb]:
                csv_writer.write("%s|%s|%s\n"%(emdb_pdb, key, map_blob_statistics[emdb_pdb][key]))

            csv_writer.close()
        
        except Exception as e:
            print("Error occured during {}".format(emdb_pdb))
            print(e)
            print(e.args)
            csv_writer.close()
        
        print("#################################################################################################### \n")
        print("#################################################################################################### \n")
    
    with open(pickle_output_file, "wb") as result_file:
        pickle.dump(map_blob_statistics, result_file)
        
def plot_linear_regression(data_input, x_col, y_col, x_label=None, y_label=None, title_text=None):
    import matplotlib.pyplot as plt
    def linear(x,a,b):
        return a * x + b
    from matplotlib.offsetbox import AnchoredText
    f, ax = plt.subplots(1,1)

    data_unsort = data_input.copy()
    data=data_unsort.sort_values(by=x_col)
    x_data = data[x_col]
    y_data = data[y_col]
    
    
    
    ax.plot(x_data, y_data,'bo')
    ax.plot(x_data, x_data, 'r-')
    equation = "y = x"
    legend_text = equation
    anchored_text=AnchoredText(legend_text, loc=2)
    ax.add_artist(anchored_text)
    if x_label is not None:
        ax.set_xlabel(x_label)
    else:
        ax.set_xlabel(x_col)
        
    if y_label is not None:
        ax.set_ylabel(y_label)
    else:
        ax.set_ylabel(y_col)
    ax.set_title(title_text)    

        
def analyse_pickle_file(result_type='mean_surface_to_volume'):
    import pickle
    import os
    import statistics
    import pandas as pd
    
    parent_folder = "/mnt/c/Users/abharadwaj1/Downloads/ForUbuntu/LocScale/Test/large_scale_quality_study"
    pickle_output_file = os.path.join(parent_folder, "surface_analysis_combined.pickle")
    
    with open(pickle_output_file, 'rb') as output_file:
        map_blob_statistics = pickle.load(output_file)
    
    
    statistics_emdb = {}
    
    
    for emd_pdb in map_blob_statistics.keys():
        unsharpened_map_blob_stats = map_blob_statistics[emd_pdb]['unsharpened_map_blob_statistics']
        sharpened_map_blob_stats = map_blob_statistics[emd_pdb]['sharpened_blob_statistics']
        
        unsharp_surface_area_to_volume_array = unsharpened_map_blob_stats['surface_area_to_volume_array']
        unsharp_surface_area_array = unsharpened_map_blob_stats['surface_area_array']
        unsharp_volume_array = unsharpened_map_blob_stats['volume_array']
        unsharp_global_surface_area_to_volume = unsharpened_map_blob_stats['total_surface_area'] / unsharpened_map_blob_stats['total_volume']

        sharp_surface_area_to_volume_array = sharpened_map_blob_stats['surface_area_to_volume_array']
        sharp_surface_area_array = sharpened_map_blob_stats['surface_area_array']
        sharp_volume_array = sharpened_map_blob_stats['volume_array']
        sharp_global_surface_area_to_volume = sharpened_map_blob_stats['total_surface_area'] / sharpened_map_blob_stats['total_volume']
        
     
        stats = {}
        ## Statistics
        # > Mean
        # > Median - 50 percentile 
        # > Maximum
        # > 
        
        ### Gather stats for unsharp 
        stats['unsharp_surface_to_volume_global'] = unsharp_global_surface_area_to_volume
        stats['unsharp_mean_surface_to_volume'] = unsharp_surface_area_to_volume_array.mean()
        stats['unsharp_median_surface_to_volume'] = statistics.median(unsharp_surface_area_to_volume_array)
        stats['unsharp_max_surface_to_volume'] = unsharp_surface_area_to_volume_array.max()
        
        stats['unsharp_mean_surface_area'] = unsharp_surface_area_array.mean()
        stats['unsharp_median_surface_area'] = statistics.median(unsharp_surface_area_array)
        stats['unsharp_max_surface_area'] = unsharp_surface_area_array.max()
        
        stats['unsharp_mean_volume'] = unsharp_volume_array.mean()
        stats['unsharp_median_volume'] = statistics.median(unsharp_volume_array)
        stats['unsharp_max_volume'] = unsharp_volume_array.max()
        
        ### Gather stats for sharpened map
        stats['sharp_surface_to_volume_global'] = sharp_global_surface_area_to_volume
        stats['sharp_mean_surface_to_volume'] = sharp_surface_area_to_volume_array.mean()
        stats['sharp_median_surface_to_volume'] = statistics.median(sharp_surface_area_to_volume_array)
        stats['sharp_max_surface_to_volume'] = sharp_surface_area_to_volume_array.max()
        
        stats['sharp_mean_surface_area'] = sharp_surface_area_array.mean()
        stats['sharp_median_surface_area'] = statistics.median(sharp_surface_area_array)
        stats['sharp_max_surface_area'] = sharp_surface_area_array.max()
        
        stats['sharp_mean_volume'] = sharp_volume_array.mean()
        stats['sharp_median_volume'] = statistics.median(sharp_volume_array)
        stats['sharp_max_volume'] = sharp_volume_array.max()        
        
        statistics_emdb[emd_pdb] = stats
        
    if result_type == 'mean_surface_to_volume':
        x = 'unsharp_mean_surface_to_volume'
        y = 'sharp_mean_surface_to_volume'
    elif result_type == 'median_surface_to_volume':
        x = 'unsharp_median_surface_to_volume'
        y = 'sharp_median_surface_to_volume'
    elif result_type == 'max_surface_to_volume':
        x = 'unsharp_max_surface_to_volume'
        y = 'sharp_max_surface_to_volume'
    elif result_type == 'mean_surface':
        x = 'unsharp_mean_surface_area'
        y = 'sharp_mean_surface_area'
    elif result_type == 'median_surface':
        x = 'unsharp_median_surface_area'
        y = 'sharp_median_surface_area'
    elif result_type == 'max_surface':
        x = 'unsharp_max_surface_area'
        y = 'sharp_max_surface_area'
    elif result_type == 'mean_volume':
        x = 'unsharp_mean_volume'
        y = 'sharp_mean_volume'
    elif result_type == 'median_volume':
        x = 'unsharp_median_volume'
        y = 'sharp_median_volume'
    elif result_type == 'max_volume':
        x = 'unsharp_max_volume'
        y = 'sharp_max_volume'
    elif result_type == "global_surface_to_volume":
        x = 'unsharp_surface_to_volume_global'
        y = 'sharp_surface_to_volume_global'
    else:
        print("unknown input, using mean_surface_to_volume as default")
        result_type = 'mean_surface_to_volume'
        x = 'unsharp_mean_surface_to_volume'
        y = 'sharp_mean_surface_to_volume'
    
    metrics_list_unique = list(list(statistics_emdb.values())[-1].keys())
    df = pd.DataFrame(data=statistics_emdb.values(), columns=list(metrics_list_unique))
    df = df.astype(float)
    

    sharpening_success = df[y] > df[x]
    sharpening_success_rate = sharpening_success.sum() / len(sharpening_success.keys()) * 100
 
    plot_linear_regression(data_input=df, x_col=x, y_col=y, title_text=result_type+" with sharpening rate {}".format(round(sharpening_success_rate,2)))
    

    print("Sharpening rate = {}% with N = {}".format(round(sharpening_success_rate,2), len(sharpening_success.keys())))
    

if __name__ == "__main__":
    
    '''
    import os
    import multiprocessing
    from locscale.utils.scaling_tools import split_sequence_evenly
    
    parent_folder = "/home/abharadwaj1/dev/new_tests_locscale/large_scale_analysis/mapdata_copy/mapdata"
    
    EMDB_PDB_ids = ["0026_6gl7", "0038_6gml", "0071_6gve", "0093_6gyn", "0094_6gyo", "0132_6h3c", "0234_6hjn", "0408_6nbd", "0415_6nbq", "4288_6fo2", "0452_6nmi", "0490_6nr8", "0492_6nra", "0567_6o0h", "0589_6nmi", "0592_6o1m", "0665_6oa9", "0776_6ku9", "10049_6rx4", "10069_6s01", "10100_6s5t", "10105_6s6t", "10106_6s6u", "10273_6sof", "10279_6sp2", "10324_6swe", "10333_6swy", "10418_6t9n", "10534_6tni", "10585_6ttu", "10595_6tut", "10617_6xt9", "20145_6oo4", "20146_6oo5", "20189_6osy", "20234_6p19", "20249_6p4h", "20254_6p5a", "20259_6p62", "20270_6p7v", "20271_6p7w", "20352_6pik", "20521_6pxm", "20986_6v0b", "21012_6v1i", "21107_6v8o", "21144_6vbu", "21391_6vv5", "3661_5no2", "3662_5no3", "3802_5of4", "3885_6el1", "3908_6eoj", "4032_5lc5", "4073_5lmn", "4074_5lmo", "4079_5lmt", "4148_5m3m", "4162_6ezo", "4192_6f6w", "4214_6fai", "4241_6fe8", "4272_6fki", "4401_6i2x", "4404_6i3m", "4429_6i84", "4588_6qm5", "4589_6qm6", "4593_6qma", "4728_6r5k", "4746_6r7x", "4759_6r8f", "4888_6ric", "4889_6rid", "4890_6rie", "4907_6rkd", "4917_6rla", "4918_6rlb", "4941_6rn3", "4983_6rqj", "7009_6ave", "7041_6b3q", "7065_6b7y", "7090_6bf6", "7334_6c23", "7335_6c24", "8911_6dt0", "8958_6e1n", "8960_6e1p", "9258_6muw", "9259_6mux", "9931_6k7g", "9934_6k7i", "9935_6k7j", "9939_6k7l", "9941_6k7m", "9695_6iok", "0193_6hcg", "0257_6hra", "0264_6hs7", "0499_6nsk", "10401_6t8h", "20449_6pqo", "20849_6uqk", "4611_6qp6", "4646_6qvb", "4733_6r69", "4789_6rb9", "7133_6bqv", "7882_6dg7", "8069_5i08", "9112_6mgv", "9298_6mzc", "9374_6nhv", "0282_6huo", "0311_6hz5", "0560_6nzu", "10365_6t23", "20220_6oxl", "20226_6p07", "3545_5mqf", "4141_5m1s", "4531_6qdw", "4571_6qk7", "4997_6rtc", "7127_6bpq", "7573_6crv", "8702_5vkq", "9610_6adq"]
    fsc_resolutions_list = ["6gl7_6.3",	"6gml_3.2",	"6gve_3.9",	"6gyn_3.4",	"6gyo_3.4",	"6h3c_3.9",	"6hjn_3.3",	"6nbd_3.2",	"6nbq_3.1",	"6fo2_4.4",	"6nmi_3.7",	"6nr8_7.8",	"6nra_7.7",	"6o0h_3.67",	"6nmi_3.7",	"6o1m_3.15",	"6oa9_3.9",	"6ku9_2.67",	"6rx4_3.3",	"6s01_3.2",	"6s5t_4.15",	"6s6t_4.1",	"6s6u_3.5",	"6sof_4.3",	"6sp2_3.33",	"6swe_3.1",	"6swy_3.2",	"6t9n_2.96",	"6tni_3.4",	"6ttu_3.7",	"6tut_3.25",	"6xt9_3.8",	"6oo4_3.3",	"6oo5_4.2",	"6osy_4.3",	"6p19_3.8",	"6p4h_3.2",	"6p5a_3.6",	"6p62_3.57",	"6p7v_4",	"6p7w_4.1",	"6pik_7.8",	"6pxm_2.1",	"6v0b_4.1",	"6v1i_3.8",	"6v8o_3.07",	"6vbu_3.1",	"6vv5_3.5",	"5no2_5.16",	"5no3_5.16",	"5of4_4.4",	"6el1_6.1",	"6eoj_3.55",	"5lc5_4.35",	"5lmn_3.55",	"5lmo_4.3",	"5lmt_4.15",	"5m3m_4",	"6ezo_4.1",	"6f6w_3.81",	"6fai_3.4",	"6fe8_4.1",	"6fki_4.3",	"6i2x_3.35",	"6i3m_3.93",	"6i84_4.4",	"6qm5_3.6",	"6qm6_3.7",	"6qma_3.7",	"6r5k_4.8",	"6r7x_3.47",	"6r8f_3.8",	"6ric_2.8",	"6rid_2.9",	"6rie_3.1",	"6rkd_3.2",	"6rla_3.9",	"6rlb_4.5",	"6rn3_4",	"6rqj_3.5",	"6ave_3.7",	"6b3q_3.7",	"6b7y_6.5",	"6bf6_6.5",	"6c23_3.9",	"6c24_3.5",	"6dt0_3.7",	"6e1n_3.7",	"6e1p_3.7",	"6muw_3.6",	"6mux_3.9",	"6k7g_3.3",	"6k7i_3.22",	"6k7j_3.08",	"6k7l_2.83",	"6k7m_2.95",	"6iok_3.64",	"6hcg_4.3",	"6hra_3.7",	"6hs7_4.6",	"6nsk_2.7",	"6t8h_3.77",	"6pqo_2.88",	"6uqk_3.77",	"6qp6_3.2",	"6qvb_4.34",	"6r69_3.65",	"6rb9_3.2",	"6bqv_3.1",	"6dg7_3.32",	"5i08_4.04",	"6mgv_3.1",	"6mzc_4.5",	"6nhv_3.5",	"6huo_3.26",	"6hz5_4.2",	"6nzu_3.2",	"6t23_3.1",	"6oxl_3.5",	"6p07_3.2",	"5mqf_5.9",	"5m1s_6.7",	"6qdw_2.83",	"6qk7_3.3",	"6rtc_3.96",	"6bpq_4.1",	"6crv_3.2","5vkq_3.55","6adq_3.5"]
        
    processes = 20
    jobs = []
    split_sequences = split_sequence_evenly(EMDB_PDB_ids, processes)
    
    for i in range(processes):
        list_of_sequence_for_this_process = split_sequences[i]
        output_pdf_folder = "/home/abharadwaj1/dev/locscale/locscale/analysis/threshold"
        output_filename = "threshold_analysis_process_{}.pdf".format(i)
        pdf_filename = os.path.join(output_pdf_folder, output_filename)
        process = multiprocessing.Process(target=get_threshold_analysis_data, args=(parent_folder, list_of_sequence_for_this_process, 
                                                                 pdf_filename))
        jobs.append(process)
    
    for job in jobs:
        job.start()
    
    for job in jobs:
        job.join()
    
    print("Threshold analysis complete!")
    '''
    
