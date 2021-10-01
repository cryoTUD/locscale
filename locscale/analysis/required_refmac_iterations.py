#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep  8 18:27:47 2021

@author: alok
"""
import gemmi
import mrcfile
import matplotlib.pyplot as plt
import numpy as np
from locscale_headers import run_refmac
from emmer.pdb.pdb_utils import get_bfactors

## Script to analyse how many REFMAC iterations are required for sufficient refinement

input_model_path = "/mnt/c/Users/abharadwaj1/Downloads/ForUbuntu/LocScale/test_locscale/emd0038/pdb6gml.ent"

input_map_path = "/mnt/c/Users/abharadwaj1/Downloads/ForUbuntu/LocScale/test_locscale/emd0038/emd_0038_half_map_1.map"

#validate_map_path = input_map_path
validate_map_path = "/mnt/c/Users/abharadwaj1/Downloads/ForUbuntu/LocScale/test_locscale/emd0038/emd_0038_half_map_2.map"

resolution = 3.2

gemmi_st = gemmi.read_structure(input_model_path)
dims = [gemmi_st.cell.a, gemmi_st.cell.b, gemmi_st.cell.c]

refined_bfactors_list = {}
refined_model_path = {}
for num_iter in np.arange(start=1, stop=20, step=2):
    print("Using num_iter of {} \n".format(num_iter))
    refined_model_path[num_iter] = run_refmac(model_path=input_model_path, model_name = input_model_path[:-4], map_path=input_map_path, resolution=resolution, maskdims=dims, num_iter=num_iter, verbose=True)
    
    gemmi.read_structure(refined_model_path[num_iter]).write_pdb(refined_model_path[num_iter][:-4]+"_iter"+str(num_iter)+".pdb")
    refined_bfactors_list[num_iter] = get_bfactors(refined_model_path[num_iter], return_as_list=True)
    
    plt.hist(refined_bfactors_list[num_iter],bins=100);
    print(np.histogram(refined_bfactors_list[num_iter],bins=20))


import seaborn as sns
for bfactorlist in refined_bfactors_list.values():
    sns.kdeplot(bfactorlist,bw=0.1)
plt.yscale("log")
plt.savefig("bfactor_histograms.jpg")

plt.figure(2)
emmap = mrcfile.open(input_map_path).data
from emmer.pdb.pdb_to_map import pdb2map
from emmer.ndimage.map_tools import compute_real_space_correlation as rsc
rscc = {}
for it in refined_model_path.keys():
    pdb = refined_model_path[it][:-4]+"_iter"+str(it)+".pdb"
    simmap = pdb2map(input_pdb=pdb, apix=1.2156, size=(256,256,256))
    rscc[it] = rsc(simmap, emmap)

plt.plot(list(rscc.keys()),list(rscc.values()))
plt.yscale("linear")
plt.savefig("RSCC.jpg")

import pickle
with open("bfactor_histograms.pickle","wb") as file:
    pickle.dump(refined_bfactors_list,file)




    
    