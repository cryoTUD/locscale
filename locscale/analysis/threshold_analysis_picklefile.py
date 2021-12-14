#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec 10 19:48:42 2021

@author: alok
"""


import mrcfile
import numpy as np
import csv
import pickle
import os
import matplotlib.pyplot as plt

folder = "/mnt/c/Users/abharadwaj1/Downloads/ForUbuntu/LocScale/tests/threshold_analysis"
pickle_filename = os.path.join(folder, "threshold_analysis_combined.pickle")

def plot_dictionary(dictionary):
    x = np.array(list(dictionary.keys()))
    y = np.array(list(dictionary.values()))
    plt.plot(x, y, 'k.-')
    
with open(pickle_filename, "rb") as file:
    threshold_analysis = pickle.load(file)
    
for emdb_pdb in threshold_analysis.keys():
    unsharp_result = threshold_analysis[emdb_pdb]['unsharp_analysis']
    sharp_result = threshold_analysis[emdb_pdb]['sharp_analysis']
    
    