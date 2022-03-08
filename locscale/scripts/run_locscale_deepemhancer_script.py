#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar  8 13:13:39 2022

@author: alok
"""

import mrcfile
import numpy as np
import csv
import pickle
import os
from subprocess import run
from locscale.utils.map_quality import map_quality_kurtosis, map_quality_pdb_multiple, measure_debye_pwlf
from locscale.include.emmer.ndimage.map_quality_tools import calculate_adjusted_surface_area, calculate_unit_surface_area
from locscale.include.emmer.pdb.pdb_tools import find_wilson_cutoff
from locscale.pseudomodel.pseudomodel_headers import number_of_segments
from locscale.utils.webscrape_symmetry import extract_symmetry
 

def run_locscale_script(parent_folder, output_folder, EMDB_PDB_ids, fsc_resolutions_list, process_id=None):
    
    
    locscale_script = os.path.join(os.environ["LOCSCALE_PATH"], "locscale","main.py")
    
    symmetry_dictionary = {
        "10279":"C6",
        "4192":"C1",
        "10106":"C2",
        "3908":"C1",
        "21107":"C1",
        "4162":"C2",
        "8960":"C2",
        "0234":"C3",
        "9610":"C2",
        "8702":"C4",
        "0193":"C15",
        "0408":"C2",
        "4272":"C1"}
    
    fsc_resolution = {}
    for data in fsc_resolutions_list:
        PDB_ID_FSC_RESOLUTION = data.split("_")[0]
        resolution = float(data.split("_")[1])
        fsc_resolution[PDB_ID_FSC_RESOLUTION] = resolution
        
    halfmap_prefix = "EMDBmaps/halfmaps"
    mask_prefix = "EMDBmaps/"
        
    for emdb_pdb in EMDB_PDB_ids:
        print("#################################################################################################### \n")
        
        print(emdb_pdb)
        try:
            emdb_id = emdb_pdb.split("_")[0]
            pdb_id = emdb_pdb.split("_")[1]
            
            
            
            ## Skip items not in the symmetry dictionary
            if str(emdb_id) not in list(symmetry_dictionary.keys()):
                continue
            
            try:
                fsc_resolution_map = fsc_resolution[pdb_id]
                print("FSC resolution = {}".format(fsc_resolution_map))
            
            except:
                print("Unable to fetch FSC resolution for {}".format(emdb_pdb))
                
            
            halfmap1_name = "emd_{}_half_map_1.map".format(emdb_id)
            halfmap2_name = "emd_{}_half_map_2.map".format(emdb_id)
            mask_filename = "emd_{}_confidence.map".format(emdb_id)
            
            
            halfmap1_path = os.path.join(parent_folder, emdb_pdb, halfmap_prefix, halfmap1_name)
            halfmap2_path = os.path.join(parent_folder, emdb_pdb, halfmap_prefix, halfmap2_name)
            mask_path = os.path.join(parent_folder, emdb_pdb, mask_prefix, mask_filename)
            
            halfmap_1_scaled_name = "locscale_{}_halfmap_1.mrc".format(emdb_id)
            halfmap_2_scaled_name = "locscale_{}_halfmap_2.mrc".format(emdb_id)
            
            halfmap_1_scaled_path = os.path.join(output_folder, halfmap_1_scaled_name)
            halfmap_2_scaled_path = os.path.join(output_folder, halfmap_2_scaled_name)
            
            processing_files_1_path = os.path.join(output_folder, "processing_files_1_{}".format(emdb_id))
            processing_files_2_path = os.path.join(output_folder, "processing_files_2_{}".format(emdb_id))
            
            output_script_1 = open(os.path.join (output_folder,"locscale_output_{}.txt".format(emdb_pdb)),"w")
            output_script_2 = open(os.path.join(output_folder,"locscale_output_{}.txt".format(emdb_pdb)),"w")
            symmetry = symmetry_dictionary[str(emdb_id)]
            
            print("Now running Locscale on Halfmap 1")
            locscale_1_command = ["mpirun","-np","4","python",locscale_script,"-em",halfmap1_path, "-ma",mask_path,"-res",str(fsc_resolution_map),"--symmetry",symmetry,"-o",
                                  halfmap_1_scaled_path,"-op",processing_files_1_path,"--verbose","--mpi"]
            
            
            print(" ".join(locscale_1_command))
            run(locscale_1_command, stdout=output_script_1)
            
            print("Now running Locscale on Halfmap 2")
            locscale_2_command = ["mpirun","-np","4","python",locscale_script,"-em",halfmap2_path, "-ma",mask_path,"-res",str(fsc_resolution_map),"--symmetry",symmetry,"-o",
                                  halfmap_2_scaled_path,"-op",processing_files_2_path,"--verbose","--mpi"]
            print(" ".join(locscale_2_command))
            run(locscale_2_command, stdout=output_script_2)
            
            
            output_script_1.close()
            output_script_2.close()
            
           
        except Exception as e:
            print("Error occured during {}".format(emdb_pdb))
            print(e)
            print(e.args)
        
            output_script_1.close()
            output_script_2.close()
            
        
        print("#################################################################################################### \n")
        print("#################################################################################################### \n")


if __name__ == "__main__":
    
    import multiprocessing
    from locscale.utils.scaling_tools import split_sequence_evenly
    import os
    
    parent_folder = "/home/abharadwaj1/dev/large_scale_analysis_locscale/mapdata"
    output_folder = "/home/abharadwaj1/dev/large_scale_analysis_locscale/halfmap_sharpening_deepEmhancer_dataset"
    
    EMDB_PDB_ids = ["0026_6gl7", "0038_6gml", "0071_6gve", "0093_6gyn", "0094_6gyo", "0132_6h3c", "0234_6hjn", "0408_6nbd", "0415_6nbq", "4288_6fo2", "0452_6nmi", "0490_6nr8", "0492_6nra", "0567_6o0h", "0589_6nmi", "0592_6o1m", "0665_6oa9", "0776_6ku9", "10049_6rx4", "10069_6s01", "10100_6s5t", "10105_6s6t", "10106_6s6u", "10273_6sof", "10279_6sp2", "10324_6swe", "10333_6swy", "10418_6t9n", "10534_6tni", "10585_6ttu", "10595_6tut", "10617_6xt9", "20145_6oo4", "20146_6oo5", "20189_6osy", "20234_6p19", "20249_6p4h", "20254_6p5a", "20259_6p62", "20270_6p7v", "20271_6p7w", "20352_6pik", "20521_6pxm", "20986_6v0b", "21012_6v1i", "21107_6v8o", "21144_6vbu", "21391_6vv5", "3661_5no2", "3662_5no3", "3802_5of4", "3885_6el1", "3908_6eoj", "4032_5lc5", "4073_5lmn", "4074_5lmo", "4079_5lmt", "4148_5m3m", "4162_6ezo", "4192_6f6w", "4214_6fai", "4241_6fe8", "4272_6fki", "4401_6i2x", "4404_6i3m", "4429_6i84", "4588_6qm5", "4589_6qm6", "4593_6qma", "4728_6r5k", "4746_6r7x", "4759_6r8f", "4888_6ric", "4889_6rid", "4890_6rie", "4907_6rkd", "4917_6rla", "4918_6rlb", "4941_6rn3", "4983_6rqj", "7009_6ave", "7041_6b3q", "7065_6b7y", "7090_6bf6", "7334_6c23", "7335_6c24", "8911_6dt0", "8958_6e1n", "8960_6e1p", "9258_6muw", "9259_6mux", "9931_6k7g", "9934_6k7i", "9935_6k7j", "9939_6k7l", "9941_6k7m", "9695_6iok", "0193_6hcg", "0257_6hra", "0264_6hs7", "0499_6nsk", "10401_6t8h", "20449_6pqo", "20849_6uqk", "4611_6qp6", "4646_6qvb", "4733_6r69", "4789_6rb9", "7133_6bqv", "7882_6dg7", "8069_5i08", "9112_6mgv", "9298_6mzc", "9374_6nhv", "0282_6huo", "0311_6hz5", "0560_6nzu", "10365_6t23", "20220_6oxl", "20226_6p07", "3545_5mqf", "4141_5m1s", "4531_6qdw", "4571_6qk7", "4997_6rtc", "7127_6bpq", "7573_6crv", "8702_5vkq", "9610_6adq"]
    fsc_resolutions_list = ["6gl7_6.3",	"6gml_3.2",	"6gve_3.9",	"6gyn_3.4",	"6gyo_3.4",	"6h3c_3.9",	"6hjn_3.3",	"6nbd_3.2",	"6nbq_3.1",	"6fo2_4.4",	"6nmi_3.7",	"6nr8_7.8",	"6nra_7.7",	"6o0h_3.67",	"6nmi_3.7",	"6o1m_3.15",	"6oa9_3.9",	"6ku9_2.67",	"6rx4_3.3",	"6s01_3.2",	"6s5t_4.15",	"6s6t_4.1",	"6s6u_3.5",	"6sof_4.3",	"6sp2_3.33",	"6swe_3.1",	"6swy_3.2",	"6t9n_2.96",	"6tni_3.4",	"6ttu_3.7",	"6tut_3.25",	"6xt9_3.8",	"6oo4_3.3",	"6oo5_4.2",	"6osy_4.3",	"6p19_3.8",	"6p4h_3.2",	"6p5a_3.6",	"6p62_3.57",	"6p7v_4",	"6p7w_4.1",	"6pik_7.8",	"6pxm_2.1",	"6v0b_4.1",	"6v1i_3.8",	"6v8o_3.07",	"6vbu_3.1",	"6vv5_3.5",	"5no2_5.16",	"5no3_5.16",	"5of4_4.4",	"6el1_6.1",	"6eoj_3.55",	"5lc5_4.35",	"5lmn_3.55",	"5lmo_4.3",	"5lmt_4.15",	"5m3m_4",	"6ezo_4.1",	"6f6w_3.81",	"6fai_3.4",	"6fe8_4.1",	"6fki_4.3",	"6i2x_3.35",	"6i3m_3.93",	"6i84_4.4",	"6qm5_3.6",	"6qm6_3.7",	"6qma_3.7",	"6r5k_4.8",	"6r7x_3.47",	"6r8f_3.8",	"6ric_2.8",	"6rid_2.9",	"6rie_3.1",	"6rkd_3.2",	"6rla_3.9",	"6rlb_4.5",	"6rn3_4",	"6rqj_3.5",	"6ave_3.7",	"6b3q_3.7",	"6b7y_6.5",	"6bf6_6.5",	"6c23_3.9",	"6c24_3.5",	"6dt0_3.7",	"6e1n_3.7",	"6e1p_3.7",	"6muw_3.6",	"6mux_3.9",	"6k7g_3.3",	"6k7i_3.22",	"6k7j_3.08",	"6k7l_2.83",	"6k7m_2.95",	"6iok_3.64",	"6hcg_4.3",	"6hra_3.7",	"6hs7_4.6",	"6nsk_2.7",	"6t8h_3.77",	"6pqo_2.88",	"6uqk_3.77",	"6qp6_3.2",	"6qvb_4.34",	"6r69_3.65",	"6rb9_3.2",	"6bqv_3.1",	"6dg7_3.32",	"5i08_4.04",	"6mgv_3.1",	"6mzc_4.5",	"6nhv_3.5",	"6huo_3.26",	"6hz5_4.2",	"6nzu_3.2",	"6t23_3.1",	"6oxl_3.5",	"6p07_3.2",	"5mqf_5.9",	"5m1s_6.7",	"6qdw_2.83",	"6qk7_3.3",	"6rtc_3.96",	"6bpq_4.1",	"6crv_3.2","5vkq_3.55","6adq_3.5"]
        
    processes = 4
    jobs = []
    split_sequences = split_sequence_evenly(EMDB_PDB_ids, processes)
    
    for i in range(processes):
        output_folder_process = os.path.join(output_folder, "process_{}".format(i))
        if not os.path.exists(output_folder_process):
            os.mkdir(output_folder_process)
        
        list_of_sequence_for_this_process = split_sequences[i]
        process_id = i
        process = multiprocessing.Process(target=run_locscale_script,
                                          args=(parent_folder, output_folder_process, list_of_sequence_for_this_process, fsc_resolutions_list, process_id))
        jobs.append(process)
    
    for job in jobs:
        job.start()
    
    for job in jobs:
        job.join()
    
    print("LocScale analysis complete!")

    
    
    