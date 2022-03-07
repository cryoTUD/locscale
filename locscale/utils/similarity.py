#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct  4 22:27:18 2021

@author: alok
"""
import numpy as np
import mrcfile


def similarity_threshold(reference_map, threshold_reference, target_map, num_bins=100, min_threshold_target=None):
    from emmer.ndimage.map_tools import compute_real_space_correlation as rsc
    from tqdm import tqdm
    
    
    
    def binarizeMap(emmap, reference):
        return (emmap>reference).astype(np.int_)
    
    if min_threshold_target is None:
        min_threshold_target = 0
    else:
        min_threshold_target = min_threshold_target
    threshold_range_target = np.linspace(min_threshold_target, target_map.max(), num_bins)
    
    binarised_reference = binarizeMap(reference_map,threshold_reference)
    rscc = {}
    for threshold in tqdm(threshold_range_target, desc="Thresholding..."):
        binarised_target = binarizeMap(target_map, threshold)
        rscc[threshold] = rsc(binarised_reference, binarised_target)
    
    return rscc

def find_threshold(reference_map_path, reference_threshold, list_of_target_map_paths, mask_path=None, num_bins=100,min_threshold_target=None):
    import matplotlib.pyplot as plt
    
    reference_map = mrcfile.open(reference_map_path).data
    if mask_path is not None:
        mask = mrcfile.open(mask_path).data
        reference_map = reference_map * mask
    rscc_curve = {}
    target_threshold = {}
    legend_name = []
    for target_map_path in list_of_target_map_paths:
        map_name = target_map_path.split("/")[-1]
        target_map = mrcfile.open(target_map_path).data
        if mask_path is not None:
            target_map = target_map * mask
        
        rsc_curve = similarity_threshold(reference_map, threshold_reference=reference_threshold, target_map=target_map, min_threshold_target=min_threshold_target)
        rscc_curve[map_name] = rsc_curve
        target_threshold[map_name] = max(rsc_curve, key=rsc_curve.get)
        plt.plot(rsc_curve.keys(), rsc_curve.values())
        legend_name.append(map_name)
    plt.legend(legend_name)
    
    print(target_threshold)
    return target_threshold, rscc_curve

    