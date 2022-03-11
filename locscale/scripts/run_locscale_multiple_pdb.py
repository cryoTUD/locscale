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

rmsd_magnitudes = [0, 1, 5, 10, 15, 20]

pdb_paths = {}
locscale_paths = {}
processing_files_folder = {}

for rmsd in rmsd_magnitudes:    
    pdb_paths[rmsd] = os.path.join(input_folder,"scattered_pdbs_using_mask","pdb6y5a_rmsd_{}_pm_perturbed_using_mask.pdb".format(rmsd))
    locscale_paths[rmsd] = os.path.join(output_folder,"locscale_rmsd_{}_A.mrc".format(rmsd))
    processing_files_folder[rmsd] = os.path.join(output_folder, "processing_files_rmsd_{}_A".format(rmsd))
    

print("output locscale maps: \n")
print(locscale_paths.values())
print("output processing files folder: \n")
print(processing_files_folder.values())

resolution = "2.8"
locscale_script = os.path.join(os.environ["LOCSCALE_PATH"], "locscale","main.py")

print("Running the following scripts: \n")

locscale_run_scripts = {}
for rmsd in rmsd_magnitudes:
    print("\n################\n")
    script = ["mpirun","-np","4","python",locscale_script,"--half_map1",halfmap_1_path, "--half_map2",halfmap_2_path,"--mask",mask_path,
              "--ref_resolution",resolution, "--verbose","--mpi",
              "--model_coordinates",pdb_paths[rmsd],"--outfile",locscale_paths[rmsd],
              "--output_processing_files",processing_files_folder[rmsd],"--skip_refine", "--output_report"]
    
    print(" ".join(script))
    locscale_run_scripts[rmsd] = script

print("\n############### START ###################\n")


start_program = True
if start_program:
    for rmsd in rmsd_magnitudes:
        return_command = run(locscale_run_scripts[rmsd])
    
        print("===========Locscale run finished for rmsd {} A ====== ".format(rmsd))
        print("Exit code: ",return_command.returncode)
    




