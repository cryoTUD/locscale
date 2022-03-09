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
mask_path = os.path.join(input_folder, "emd_10692_confidenceMap.mrc")

perturb_magnitudes = [1, 5, 10, 15, 20]

pdb_paths = {}
locscale_paths = {}
processing_files_folder = {}

for perturb in perturb_magnitudes:    
    pdb_paths[perturb] = os.path.join(input_folder,"scattered_pdbs","pdb6y5a_rmsd_{}_A.pdb".format(perturb))
    locscale_paths[perturb] = os.path.join(output_folder,"locscale_rmsd_{}_A.mrc".format(perturb))
    processing_files_folder[perturb] = os.path.join(output_folder, "processing_files_rmsd_{}_A".format(perturb))
    

print("output locscale maps: \n")
print(locscale_paths.values())
print("output processing files folder: \n")
print(processing_files_folder.values())

resolution = "2.8"
locscale_script = os.path.join(os.environ["LOCSCALE_PATH"], "locscale","main.py")

print("Running the following scripts: \n")

locscale_run_scripts = {}
for perturb in perturb_magnitudes:
    print("\n################\n")
    script = ["mpirun","-np","4","python",locscale_script,"--half_map1",halfmap_1_path, "--half_map2",halfmap_2_path,"--mask",mask_path,
              "--ref_resolution",resolution, "--verbose","--mpi",
              "--model_coordinates",pdb_paths[perturb],"--outfile",locscale_paths[perturb],
              "--output_processing_files",processing_files_folder[perturb],"--skip_refine", "--output_report"]
    
    print(" ".join(script))
    locscale_run_scripts[perturb] = script

print("\n############### START ###################\n")


start_program = False
if start_program:
    for perturb in perturb_magnitudes:
        return_command = run(locscale_run_scripts[perturb])
    
        print("===========Locscale run finished for perturb {} A ====== ".format(perturb))
        print("Exit code: ",return_command.returncode)
    




