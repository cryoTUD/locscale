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


def get_data(parent_folder, EMDB_PDB_ids, fsc_resolutions_list,output_filename, process_id=None):
    
    from locscale.utils.map_quality import local_histogram_analysis
    
  
    pickle_output_file = os.path.join(parent_folder, "{}.pickle".format(output_filename))
    csv_output_file = os.path.join(parent_folder, "{}_csv.csv".format(output_filename))
    #EMDB_PDB_ids = ["0026_6gl7", "0038_6gml", "0071_6gve", "0093_6gyn", "0094_6gyo", "0132_6h3c", "0234_6hjn", "0408_6nbd", "0415_6nbq", "4288_6fo2", "0452_6nmi", "0490_6nr8", "0492_6nra", "0567_6o0h", "0589_6nmi", "0592_6o1m", "0665_6oa9", "0776_6ku9", "10049_6rx4", "10069_6s01", "10100_6s5t", "10105_6s6t", "10106_6s6u", "10273_6sof", "10279_6sp2", "10324_6swe", "10333_6swy", "10418_6t9n", "10534_6tni", "10585_6ttu", "10595_6tut", "10617_6xt9", "20145_6oo4", "20146_6oo5", "20189_6osy", "20234_6p19", "20249_6p4h", "20254_6p5a", "20259_6p62", "20270_6p7v", "20271_6p7w", "20352_6pik", "20521_6pxm", "20986_6v0b", "21012_6v1i", "21107_6v8o", "21144_6vbu", "21391_6vv5", "3661_5no2", "3662_5no3", "3802_5of4", "3885_6el1", "3908_6eoj", "4032_5lc5", "4073_5lmn", "4074_5lmo", "4079_5lmt", "4148_5m3m", "4162_6ezo", "4192_6f6w", "4214_6fai", "4241_6fe8", "4272_6fki", "4401_6i2x", "4404_6i3m", "4429_6i84", "4588_6qm5", "4589_6qm6", "4593_6qma", "4728_6r5k", "4746_6r7x", "4759_6r8f", "4888_6ric", "4889_6rid", "4890_6rie", "4907_6rkd", "4917_6rla", "4918_6rlb", "4941_6rn3", "4983_6rqj", "7009_6ave", "7041_6b3q", "7065_6b7y", "7090_6bf6", "7334_6c23", "7335_6c24", "8911_6dt0", "8958_6e1n", "8960_6e1p", "9258_6muw", "9259_6mux", "9931_6k7g", "9934_6k7i", "9935_6k7j", "9939_6k7l", "9941_6k7m", "9695_6iok", "0193_6hcg", "0257_6hra", "0264_6hs7", "0499_6nsk", "10401_6t8h", "20449_6pqo", "20849_6uqk", "4611_6qp6", "4646_6qvb", "4733_6r69", "4789_6rb9", "7133_6bqv", "7882_6dg7", "8069_5i08", "9112_6mgv", "9298_6mzc", "9374_6nhv", "0282_6huo", "0311_6hz5", "0560_6nzu", "10365_6t23", "20220_6oxl", "20226_6p07", "3545_5mqf", "4141_5m1s", "4531_6qdw", "4571_6qk7", "4997_6rtc", "7127_6bpq", "7573_6crv", "8702_5vkq", "9610_6adq"]
    #fsc_resolutions_list = ["6gl7_6.3",	"6gml_3.2",	"6gve_3.9",	"6gyn_3.4",	"6gyo_3.4",	"6h3c_3.9",	"6hjn_3.3",	"6nbd_3.2",	"6nbq_3.1",	"6fo2_4.4",	"6nmi_3.7",	"6nr8_7.8",	"6nra_7.7",	"6o0h_3.67",	"6nmi_3.7",	"6o1m_3.15",	"6oa9_3.9",	"6ku9_2.67",	"6rx4_3.3",	"6s01_3.2",	"6s5t_4.15",	"6s6t_4.1",	"6s6u_3.5",	"6sof_4.3",	"6sp2_3.33",	"6swe_3.1",	"6swy_3.2",	"6t9n_2.96",	"6tni_3.4",	"6ttu_3.7",	"6tut_3.25",	"6xt9_3.8",	"6oo4_3.3",	"6oo5_4.2",	"6osy_4.3",	"6p19_3.8",	"6p4h_3.2",	"6p5a_3.6",	"6p62_3.57",	"6p7v_4",	"6p7w_4.1",	"6pik_7.8",	"6pxm_2.1",	"6v0b_4.1",	"6v1i_3.8",	"6v8o_3.07",	"6vbu_3.1",	"6vv5_3.5",	"5no2_5.16",	"5no3_5.16",	"5of4_4.4",	"6el1_6.1",	"6eoj_3.55",	"5lc5_4.35",	"5lmn_3.55",	"5lmo_4.3",	"5lmt_4.15",	"5m3m_4",	"6ezo_4.1",	"6f6w_3.81",	"6fai_3.4",	"6fe8_4.1",	"6fki_4.3",	"6i2x_3.35",	"6i3m_3.93",	"6i84_4.4",	"6qm5_3.6",	"6qm6_3.7",	"6qma_3.7",	"6r5k_4.8",	"6r7x_3.47",	"6r8f_3.8",	"6ric_2.8",	"6rid_2.9",	"6rie_3.1",	"6rkd_3.2",	"6rla_3.9",	"6rlb_4.5",	"6rn3_4",	"6rqj_3.5",	"6ave_3.7",	"6b3q_3.7",	"6b7y_6.5",	"6bf6_6.5",	"6c23_3.9",	"6c24_3.5",	"6dt0_3.7",	"6e1n_3.7",	"6e1p_3.7",	"6muw_3.6",	"6mux_3.9",	"6k7g_3.3",	"6k7i_3.22",	"6k7j_3.08",	"6k7l_2.83",	"6k7m_2.95",	"6iok_3.64",	"6hcg_4.3",	"6hra_3.7",	"6hs7_4.6",	"6nsk_2.7",	"6t8h_3.77",	"6pqo_2.88",	"6uqk_3.77",	"6qp6_3.2",	"6qvb_4.34",	"6r69_3.65",	"6rb9_3.2",	"6bqv_3.1",	"6dg7_3.32",	"5i08_4.04",	"6mgv_3.1",	"6mzc_4.5",	"6nhv_3.5",	"6huo_3.26",	"6hz5_4.2",	"6nzu_3.2",	"6t23_3.1",	"6oxl_3.5",	"6p07_3.2",	"5mqf_5.9",	"5m1s_6.7",	"6qdw_2.83",	"6qk7_3.3",	"6rtc_3.96",	"6bpq_4.1",	"6crv_3.2","5vkq_3.55","6adq_3.5"]
    
    fsc_resolution = {}
    for data in fsc_resolutions_list:
        pdb_id = data.split("_")[0]
        resolution = float(data.split("_")[1])
        fsc_resolution[pdb_id] = resolution
        
    unsharpened_map_path_prefix = "EMDBmaps/"
    sharpened_map_path_prefix = "EMDBmaps/sharpened_map/"
    pdb_path_prefix = "PDBgemmis/"
    
    map_histogram_analysis = {}
    
    
    
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
            fsc_resolution_map = fsc_resolution[pdb_id]
            
            unsharpened_map_path = os.path.join(parent_folder, emdb_pdb, unsharpened_map_path_prefix, emdb_filename)
            mask_path = os.path.join(parent_folder, emdb_pdb, unsharpened_map_path_prefix, mask_filename)
            md_locscale_path = os.path.join(parent_folder, emdb_pdb, sharpened_map_path_prefix, MD_locscale_filename)
            
            
            _, unsharpened_skew_kurtosis_r2, unsharpened_mean_variance_r2 = local_histogram_analysis(unsharpened_map_path, mask_path, fsc_resolution_map)
            _, sharpened_skew_kurtosis_r2, sharpened_mean_variance_r2 = local_histogram_analysis(unsharpened_map_path, mask_path, fsc_resolution_map)
            
            
            map_histogram_analysis[emdb_pdb] = {
                'unsharpened_skew_kurtosis_r2':unsharpened_skew_kurtosis_r2, 
                'unsharpened_mean_variance_r2':unsharpened_mean_variance_r2, 
                'sharpened_skew_kurtosis_r2':sharpened_skew_kurtosis_r2, 
                'sharpened_mean_variance_r2':sharpened_mean_variance_r2, 
                }
            
            for key in map_histogram_analysis[emdb_pdb]:
                csv_writer.write("%s|%s|%s\n"%(emdb_pdb, key, map_histogram_analysis[emdb_pdb][key]))
        
            csv_writer.close()
        except Exception as e:
            print("Error occured during {}".format(emdb_pdb))
            print(e)
            print(e.args)
            csv_writer.close()
        
        print("#################################################################################################### \n")
        print("#################################################################################################### \n")
    
    with open(pickle_output_file, "wb") as result_file:
        pickle.dump(map_histogram_analysis, result_file)
        
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

        
def analyse_pickle_file(result_type='skew_kurtosis'):
    import pickle
    import os
    import statistics
    import pandas as pd
    
    parent_folder = "/mnt/c/Users/abharadwaj1/Downloads/ForUbuntu/LocScale/tests/test_quality_metrics_large"
    pickle_output_file = os.path.join(parent_folder, "histogram_analysis.pickle")
    
    with open(pickle_output_file, 'rb') as output_file:
        map_histogram_analysis = pickle.load(output_file)
    
    
    histogram_analysis_full = {}
    
    
    for emd_pdb in map_histogram_analysis.keys():
        emdb_histogram_analysis = map_histogram_analysis[emd_pdb]['unsharpened_map_histogram_analysis']

        
     
        stats = {}
        ## Statistics
        # > Mean
        # > Median - 50 percentile 
        # > Maximum
        # > 
        
        ### Gather stats for unsharp 
        
        stats['unsharpened_skew_kurtosis_r2'] = emdb_histogram_analysis['unsharpened_skew_kurtosis_r2']
        stats['unsharpened_mean_variance_r2'] = emdb_histogram_analysis['unsharpened_mean_variance_r2']
        stats['sharpened_skew_kurtosis_r2'] = emdb_histogram_analysis['sharpened_skew_kurtosis_r2']
        stats['sharpened_mean_variance_r2'] = emdb_histogram_analysis['sharpened_mean_variance_r2']
        
        histogram_analysis_full[emd_pdb] = stats
        
    if result_type == 'skew_kurtosis':
        x = 'unsharpened_skew_kurtosis_r2'
        y = 'sharpened_skew_kurtosis_r2'
    elif result_type == 'mean_variance':
        x = 'unsharpened_mean_variance_r2'
        y = 'sharpened_mean_variance_r2'
    
    else:
        print("unknown input, using skew_kurtosis as default")
        result_type = 'skew_kurtosis'
        x = 'unsharpened_skew_kurtosis_r2'
        y = 'sharpened_skew_kurtosis_r2'
    
    metrics_list_unique = list(list(histogram_analysis_full.values())[-1].keys())
    df = pd.DataFrame(data=histogram_analysis_full.values(), columns=list(metrics_list_unique))
    df = df.astype(float)
    

    sharpening_success = df[y] > df[x]
    sharpening_success_rate = sharpening_success.sum() / len(sharpening_success.keys()) * 100
 
    plot_linear_regression(data_input=df, x_col=x, y_col=y, title_text=result_type+" with sharpening rate {}".format(round(sharpening_success_rate,2)))
    

    print("Sharpening rate = {}% with N = {}".format(round(sharpening_success_rate,2), len(sharpening_success.keys())))
    

if __name__ == "__main__":
    import multiprocessing
    from locscale.utils.scaling_tools import split_sequence_evenly
    
    parent_folder = "/home/abharadwaj1/dev/new_tests_locscale/large_scale_analysis/mapdata_copy/mapdata"
    
    EMDB_PDB_ids = ["0026_6gl7", "0038_6gml", "0071_6gve", "0093_6gyn", "0094_6gyo", "0132_6h3c", "0234_6hjn", "0408_6nbd", "0415_6nbq", "4288_6fo2", "0452_6nmi", "0490_6nr8", "0492_6nra", "0567_6o0h", "0589_6nmi", "0592_6o1m", "0665_6oa9", "0776_6ku9", "10049_6rx4", "10069_6s01", "10100_6s5t", "10105_6s6t", "10106_6s6u", "10273_6sof", "10279_6sp2", "10324_6swe", "10333_6swy", "10418_6t9n", "10534_6tni", "10585_6ttu", "10595_6tut", "10617_6xt9", "20145_6oo4", "20146_6oo5", "20189_6osy", "20234_6p19", "20249_6p4h", "20254_6p5a", "20259_6p62", "20270_6p7v", "20271_6p7w", "20352_6pik", "20521_6pxm", "20986_6v0b", "21012_6v1i", "21107_6v8o", "21144_6vbu", "21391_6vv5", "3661_5no2", "3662_5no3", "3802_5of4", "3885_6el1", "3908_6eoj", "4032_5lc5", "4073_5lmn", "4074_5lmo", "4079_5lmt", "4148_5m3m", "4162_6ezo", "4192_6f6w", "4214_6fai", "4241_6fe8", "4272_6fki", "4401_6i2x", "4404_6i3m", "4429_6i84", "4588_6qm5", "4589_6qm6", "4593_6qma", "4728_6r5k", "4746_6r7x", "4759_6r8f", "4888_6ric", "4889_6rid", "4890_6rie", "4907_6rkd", "4917_6rla", "4918_6rlb", "4941_6rn3", "4983_6rqj", "7009_6ave", "7041_6b3q", "7065_6b7y", "7090_6bf6", "7334_6c23", "7335_6c24", "8911_6dt0", "8958_6e1n", "8960_6e1p", "9258_6muw", "9259_6mux", "9931_6k7g", "9934_6k7i", "9935_6k7j", "9939_6k7l", "9941_6k7m", "9695_6iok", "0193_6hcg", "0257_6hra", "0264_6hs7", "0499_6nsk", "10401_6t8h", "20449_6pqo", "20849_6uqk", "4611_6qp6", "4646_6qvb", "4733_6r69", "4789_6rb9", "7133_6bqv", "7882_6dg7", "8069_5i08", "9112_6mgv", "9298_6mzc", "9374_6nhv", "0282_6huo", "0311_6hz5", "0560_6nzu", "10365_6t23", "20220_6oxl", "20226_6p07", "3545_5mqf", "4141_5m1s", "4531_6qdw", "4571_6qk7", "4997_6rtc", "7127_6bpq", "7573_6crv", "8702_5vkq", "9610_6adq"]
    fsc_resolutions_list = ["6gl7_6.3",	"6gml_3.2",	"6gve_3.9",	"6gyn_3.4",	"6gyo_3.4",	"6h3c_3.9",	"6hjn_3.3",	"6nbd_3.2",	"6nbq_3.1",	"6fo2_4.4",	"6nmi_3.7",	"6nr8_7.8",	"6nra_7.7",	"6o0h_3.67",	"6nmi_3.7",	"6o1m_3.15",	"6oa9_3.9",	"6ku9_2.67",	"6rx4_3.3",	"6s01_3.2",	"6s5t_4.15",	"6s6t_4.1",	"6s6u_3.5",	"6sof_4.3",	"6sp2_3.33",	"6swe_3.1",	"6swy_3.2",	"6t9n_2.96",	"6tni_3.4",	"6ttu_3.7",	"6tut_3.25",	"6xt9_3.8",	"6oo4_3.3",	"6oo5_4.2",	"6osy_4.3",	"6p19_3.8",	"6p4h_3.2",	"6p5a_3.6",	"6p62_3.57",	"6p7v_4",	"6p7w_4.1",	"6pik_7.8",	"6pxm_2.1",	"6v0b_4.1",	"6v1i_3.8",	"6v8o_3.07",	"6vbu_3.1",	"6vv5_3.5",	"5no2_5.16",	"5no3_5.16",	"5of4_4.4",	"6el1_6.1",	"6eoj_3.55",	"5lc5_4.35",	"5lmn_3.55",	"5lmo_4.3",	"5lmt_4.15",	"5m3m_4",	"6ezo_4.1",	"6f6w_3.81",	"6fai_3.4",	"6fe8_4.1",	"6fki_4.3",	"6i2x_3.35",	"6i3m_3.93",	"6i84_4.4",	"6qm5_3.6",	"6qm6_3.7",	"6qma_3.7",	"6r5k_4.8",	"6r7x_3.47",	"6r8f_3.8",	"6ric_2.8",	"6rid_2.9",	"6rie_3.1",	"6rkd_3.2",	"6rla_3.9",	"6rlb_4.5",	"6rn3_4",	"6rqj_3.5",	"6ave_3.7",	"6b3q_3.7",	"6b7y_6.5",	"6bf6_6.5",	"6c23_3.9",	"6c24_3.5",	"6dt0_3.7",	"6e1n_3.7",	"6e1p_3.7",	"6muw_3.6",	"6mux_3.9",	"6k7g_3.3",	"6k7i_3.22",	"6k7j_3.08",	"6k7l_2.83",	"6k7m_2.95",	"6iok_3.64",	"6hcg_4.3",	"6hra_3.7",	"6hs7_4.6",	"6nsk_2.7",	"6t8h_3.77",	"6pqo_2.88",	"6uqk_3.77",	"6qp6_3.2",	"6qvb_4.34",	"6r69_3.65",	"6rb9_3.2",	"6bqv_3.1",	"6dg7_3.32",	"5i08_4.04",	"6mgv_3.1",	"6mzc_4.5",	"6nhv_3.5",	"6huo_3.26",	"6hz5_4.2",	"6nzu_3.2",	"6t23_3.1",	"6oxl_3.5",	"6p07_3.2",	"5mqf_5.9",	"5m1s_6.7",	"6qdw_2.83",	"6qk7_3.3",	"6rtc_3.96",	"6bpq_4.1",	"6crv_3.2","5vkq_3.55","6adq_3.5"]
        
    processes = 4
    jobs = []
    split_sequences = split_sequence_evenly(EMDB_PDB_ids, processes)
    
    for i in range(processes):
        list_of_sequence_for_this_process = split_sequences[i]
        output_filename = "histogram_analysis_process_{}".format(i)
        process = multiprocessing.Process(target=get_data, args=(parent_folder, list_of_sequence_for_this_process, 
                                                                 fsc_resolutions_list, output_filename))
        jobs.append(process)
    
    for job in jobs:
        job.start()
    
    for job in jobs:
        job.join()
    
    print("Histogram analysis complete!")


