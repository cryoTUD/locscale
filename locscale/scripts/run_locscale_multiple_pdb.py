#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar  1 20:35:01 2022

@author: alok
"""
'''
This script performs traditional locscale analysis using several PDB models as input
'''

import os
from subprocess import Popen, run,  PIPE

## Input folder contains unsharpened half maps and all perturbed PDB

input_folder = "/home/abharadwaj1/dev/faraday_tests/input_folder"
output_folder = "/home/abharadwaj1/dev/faraday_tests/output_folder"

halfmap_1_path = os.path.join(input_folder, "emd_10692_half_map_1.map")
halfmap_2_path = os.path.join(input_folder, "emd_10692_half_map_2.map")

perturb_magnitudes = [0, 1, 2, 5, 7, 10, 12, 15, 17, 20]

pdb_paths = {}
locscale_paths = {}
processing_files_folder = {}

for perturb in perturb_magnitudes:    
    pdb_paths[perturb] = os.path.join(input_folder,"pdb6y5a_perturbed_{}_A.pdb".format(perturb))
    locscale_paths[perturb] = os.path.join(output_folder,"locscale_perturbed_pdb_{}A".format(perturb))
    processing_files_folder[perturb] = os.path.join(output_folder, "processing_files_{}A".format(perturb))
    

print("output locscale maps: \n")
print(locscale_paths.values())
print("output processing files folder: \n")
print(processing_files_folder.values())

resolution = "2.8"
model_resolution = "1.8"
locscale_script = os.path.join(os.environ["LOCSCALE_PATH"], "locscale","main.py")

print("Running the following scripts: \n")

locscale_run_scripts = {}
for perturb in perturb_magnitudes:
    print("\n################\n")
    script = ["mpirun","-np","4","python",locscale_script,"--half_map1",halfmap_1_path, "--half_map2",halfmap_2_path,
              "--ref_resolution",resolution, "--model_resolution",model_resolution,"--verbose","--mpi",
              "--model_coordinates",pdb_paths[perturb],"--outfile",locscale_paths[perturb],
              "--output_processing_files",processing_files_folder[perturb]]
    
    print(" ".join(script))
    locscale_run_scripts[perturb] = script

print("\n############### START ###################\n")

for perturb in perturb_magnitudes:
    return_command = run(locscale_run_scripts[perturb], check=True)
    
    print("===========Locscale run finished for perturb {} A ====== ".format(perturb))
    print("Exit code: ",return_command.returncode)
    




