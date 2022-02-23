#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 24 13:55:06 2021

@author: alok
"""
import mrcfile
import os

def get_data():
    from locscale.include.emmer.ndimage.map_utils import compute_FDR_confidenceMap_easy
    from locscale.include.emmer.ndimage.emdb_class import EMDB
    from locscale.utils.webscrape_symmetry import extract_threshold

#    parent_folder = "/home/abharadwaj1/dev/new_tests_locscale/large_scale_analysis/mapdata_copy/mapdata"
   # pickle_output_file = os.path.join(parent_folder, "map_quality_result.pickle")
   # csv_output_file = os.path.join(parent_folder, "quality_results_csv.csv")
    EMDB_PDB_ids = ["0026_6gl7", "0038_6gml", "0071_6gve", "0093_6gyn", "0094_6gyo", "0132_6h3c", "0234_6hjn", "0408_6nbd", "0415_6nbq", "4288_6fo2", "0452_6nmi", "0490_6nr8", "0492_6nra", "0567_6o0h", "0589_6nmi", "0592_6o1m", "0665_6oa9", "0776_6ku9", "10049_6rx4", "10069_6s01", "10100_6s5t", "10105_6s6t", "10106_6s6u", "10273_6sof", "10279_6sp2", "10324_6swe", "10333_6swy", "10418_6t9n", "10534_6tni", "10585_6ttu", "10595_6tut", "10617_6xt9", "20145_6oo4", "20146_6oo5", "20189_6osy", "20234_6p19", "20249_6p4h", "20254_6p5a", "20259_6p62", "20270_6p7v", "20271_6p7w", "20352_6pik", "20521_6pxm", "20986_6v0b", "21012_6v1i", "21107_6v8o", "21144_6vbu", "21391_6vv5", "3661_5no2", "3662_5no3", "3802_5of4", "3885_6el1", "3908_6eoj", "4032_5lc5", "4073_5lmn", "4074_5lmo", "4079_5lmt", "4148_5m3m", "4162_6ezo", "4192_6f6w", "4214_6fai", "4241_6fe8", "4272_6fki", "4401_6i2x", "4404_6i3m", "4429_6i84", "4588_6qm5", "4589_6qm6", "4593_6qma", "4728_6r5k", "4746_6r7x", "4759_6r8f", "4888_6ric", "4889_6rid", "4890_6rie", "4907_6rkd", "4917_6rla", "4918_6rlb", "4941_6rn3", "4983_6rqj", "7009_6ave", "7041_6b3q", "7065_6b7y", "7090_6bf6", "7334_6c23", "7335_6c24", "8911_6dt0", "8958_6e1n", "8960_6e1p", "9258_6muw", "9259_6mux", "9931_6k7g", "9934_6k7i", "9935_6k7j", "9939_6k7l", "9941_6k7m", "9695_6iok", "0193_6hcg", "0257_6hra", "0264_6hs7", "0499_6nsk", "10401_6t8h", "20449_6pqo", "20849_6uqk", "4611_6qp6", "4646_6qvb", "4733_6r69", "4789_6rb9", "7133_6bqv", "7882_6dg7", "8069_5i08", "9112_6mgv", "9298_6mzc", "9374_6nhv", "0282_6huo", "0311_6hz5", "0560_6nzu", "10365_6t23", "20220_6oxl", "20226_6p07", "3545_5mqf", "4141_5m1s", "4531_6qdw", "4571_6qk7", "4997_6rtc", "7127_6bpq", "7573_6crv", "8702_5vkq", "9610_6adq"]
    main_folder = "/home/alok/dev/threshold_analysis"
    
    csv_output_file = os.path.join(main_folder, "fdr_threshold.csv")
    csv_writer = open(csv_output_file, mode="w")
    csv_writer.write("EMDB_ID,FDR threshold,Reporter threshold\n")
    csv_writer.close()
    
    for emdb_pdb in EMDB_PDB_ids:
        print("#################################################################################################### \n")
        csv_writer = open(csv_output_file, mode="a")
        print(emdb_pdb)
        emdb_id = emdb_pdb.split("_")[0]
        pdb_id = emdb_pdb.split("_")[1]
        DOWNLOAD_STEP=False
        FDR_STEP=False
        WEBSCRAPE_STEP=False
        try:
            emdb_filename = os.path.join(main_folder,"emd_{}.map".format(emdb_id))
            EMD = EMDB(emdb_id)
            EMD.download(main_folder)
            deposited_map = mrcfile.open(emdb_filename).data
            apix = mrcfile.open(emdb_filename).voxel_size.tolist()[0]
            length = deposited_map.shape[0]
            DOWNLOAD_STEP=True
            fdr_map, fdr = compute_FDR_confidenceMap_easy(deposited_map,apix,41, fdr=0.01, folder=main_folder)
            print("FDR - {}: {}".format(emdb_id, fdr))
            FDR_STEP=True
            reported_threshold = extract_threshold(emdb_id)
            print("Reported - {}: {}".format(emdb_id, reported_threshold))
            WEBSCRAPE_STEP = True
            
            csv_writer.write("%s,%s,%s\n"%(emdb_id, fdr, reported_threshold))
            csv_writer.close()
            if os.path.exists(emdb_filename):
                os.remove(emdb_filename)
                
        except:
            print("Problem with {}. Download = {}, FDR = {}, WEBSCRAPE = {}".format(emdb_id, DOWNLOAD_STEP, FDR_STEP, WEBSCRAPE_STEP))
            csv_writer.write("%s,%s,%s,%s\n"%(emdb_id,DOWNLOAD_STEP,FDR_STEP,WEBSCRAPE_STEP))
            csv_writer.close()
            if os.path.exists(emdb_filename):
                os.remove(emdb_filename)
                
            if not DOWNLOAD_STEP  and not FDR_STEP and not WEBSCRAPE_STEP:
                raise

    
        print("#################################################################################################### \n")

get_data()