#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Oct 16 22:16:33 2021

@author: alok
"""

## Analysis tool for Locscale parameters

import mrcfile
import os
import numpy as np
import shutil
import time
from datetime import datetime
from subprocess import run, PIPE

locscale_path = os.environ['LOCSCALE_PATH']
run_script = os.path.join(locscale_path, "locscale","main.py")

output_folder = "/home/abharadwaj1/dev/new_tests_locscale/parametric_analysis/"

input_map = os.path.join(locscale_path, "tests","test_data","emd5778_map.mrc")
input_mask = os.path.join(locscale_path, "tests","test_data","emd5778_mask.mrc")
apix = float(mrcfile.open(input_map).voxel_size.x)
input_map_res = 3.4

analysis_parameter_list = ['pm','pm_it','ref_it','add_blur','s','fdr_f','wn','mres']

parameter_values_range = {}
parameter_values_range['pseudomodel_method'] = ['gradient','kick']
parameter_values_range['total_iterations'] = list(np.arange(1,26,1,dtype=int))
parameter_values_range['distance'] = [1, 1.2, 1.4, 1.6, 1.8, 2]
parameter_values_range['refmac_iterations'] = list(np.arange(1,16,1,dtype=int))
parameter_values_range['add_blur'] = list(np.arange(5,30,5,dtype=int))
parameter_values_range['smooth_factor'] = list(np.arange(0.1,1.1,0.2).round(2))
parameter_values_range['fdr_filter'] = list(np.arange(2,10,1))
parameter_values_range['window_size'] = list(np.arange(15,45,5,dtype=int))
parameter_values_range['model_resolution'] = [3.4,round(apix*2,2)]

print("Analysis parameters:\n")
print(parameter_values_range)

os.chdir(output_folder)

for parameter in parameter_values_range.keys():
    parameter_folder_path = os.path.join(output_folder,parameter)
    os.mkdir(parameter_folder_path)
    os.chdir(parameter_folder_path)
    for value in parameter_values_range[parameter]:
        if isinstance(value, float):
            value_folder_name = str(round(value*100,2))
        else:
            value_folder_name = str(value)
        value_folder_path = os.path.join(parameter_folder_path, value_folder_name)
        
        os.mkdir(value_folder_path)
        os.chdir(value_folder_path)
        outputfile = open(os.path.join(value_folder_path,"log.txt"),"w+")
        input_map_copy = os.path.join(value_folder_path,"input_map.mrc")
        input_mask_copy = os.path.join(value_folder_path,"input_mask.mrc")
        shutil.copyfile(input_map, input_map_copy)
        shutil.copyfile(input_mask, input_mask_copy)
        locscale_run_cmd = ['mpirun','-np','4','python',run_script,'-em',input_map_copy,'-ma',input_mask_copy,'-res',str(input_map_res),'-v','-mpi']
        locscale_run_cmd.append("--"+parameter)
        locscale_run_cmd.append(str(value))
        print(" ".join(locscale_run_cmd))
        run(locscale_run_cmd, stdout=outputfile)
        
        time.sleep(1)
        print(str(datetime.now()))
        outputfile.close()
                

        

