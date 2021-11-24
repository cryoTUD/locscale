#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 26 13:46:15 2021

@author: alok
"""
## Script to calculate adjusted surface area of a map

import os
import numpy as np

def mesh_surface_area(data, threshold, apix):
    from skimage import measure
    print(data.min(), data.max())
    verts, faces,_,_ = measure.marching_cubes(data, threshold)
    surface_area = measure.mesh_surface_area(verts, faces) * apix**2
    return surface_area

def calculate_surface_area(binary_data, spacing, origin):
    from tvtk.api import tvtk
    from tvtk.common import configure_input

    grid = tvtk.ImageData(spacing=spacing, origin=origin)
    grid.point_data.scalars = binary_data.T.ravel()
    grid.point_data.scalars.name = 'scalars'
    grid.dimensions = binary_data.shape

    iso = tvtk.ImageMarchingCubes()
    configure_input(iso, grid)  
    
    mass = tvtk.MassProperties()
    configure_input(mass, iso)
    surface_area = mass.surface_area
    
    return surface_area

def calculate_surface_area_at_threshold(emmap, apix_tuple, origin, reference_threshold):
    binarised_emmap = (emmap>reference_threshold).astype(np.int_)
    surface_area = mesh_surface_area(binarised_emmap, reference_threshold, apix_tuple[0])
    return surface_area


def count_distinct_regions(emmap, reference_threshold):
    from skimage import measure
    
    binarised_emmap = (emmap>reference_threshold).astype(np.int_)
    labels, num_regions = measure.label(binarised_emmap, background=0, return_num=True)
    
    return num_regions

def find_volume_matching_threshold(emmap, reference_volume, apix, num_bins=100):
    threshold_bins = np.linspace(0, emmap.max(), num=num_bins)
    
    for threshold in threshold_bins:
        binarised_map = (emmap>=threshold).astype(np.int_)
        sum_of_voxels = binarised_map.sum()
        volume_real_units = sum_of_voxels * apix**3
        if volume_real_units <= reference_volume:
            matching_threshold = threshold
            break
    
    
    return matching_threshold

def find_c_scale(sharpest_map, most_blurred_map, reference_threshold, apix_tuple, origin):
    print("Calculating scale factor...")
    SA_blurred = calculate_surface_area_at_threshold(
        most_blurred_map, reference_threshold=reference_threshold, apix_tuple=apix_tuple, origin=origin)
    
    SA_sharpened = calculate_surface_area_at_threshold(
        sharpest_map, reference_threshold=reference_threshold, apix_tuple=apix_tuple, origin=origin)
    
    num_regions_blurred = count_distinct_regions(most_blurred_map, reference_threshold)
    num_regions_sharpened = count_distinct_regions(sharpest_map, reference_threshold)
    
    C_scale = (SA_blurred - SA_sharpened) / (num_regions_blurred - num_regions_sharpened)
    
    print("Scaling factor is: {:.2f}".format(C_scale))
    
    asa_blurred = SA_blurred - C_scale * num_regions_blurred
    asa_sharpened = SA_sharpened - C_scale * num_regions_sharpened
    

    return C_scale, asa_blurred, asa_sharpened

def apply_b_sharpen(emmap, apix, b_sharpen, d_cut, b_blur=200, k=10):
    from locscale.include.emmer.ndimage.profile_tools import compute_radial_profile, frequency_array
    from locscale.include.emmer.ndimage.map_tools import set_radial_profile, compute_scale_factors
    
    emmap_profile, radii = compute_radial_profile(emmap,return_indices=True)
    freq = frequency_array(amplitudes=emmap_profile, apix=apix)
    
    A_sharpen = np.exp(0.25 * b_sharpen * freq**2)
    A_blur = np.exp(-0.25 * b_blur * freq**2)
    
    w_sharpen = 1 / (1+np.exp(k * (d_cut - 1/freq)))
    w_blur = 1 - w_sharpen
    
    sharpening_profile = w_sharpen * A_sharpen + w_blur * A_blur
    sharpened_profile = emmap_profile * sharpening_profile

    scale_factors = compute_scale_factors(emmap_profile, sharpened_profile)
    sharpened_map = set_radial_profile(emmap, scale_factors, radii)
    
    return sharpened_map

def calculate_unit_surface_area(emmap_path, mask_path, mask_emmap=False):
    import mrcfile
    
    print("Calculating unit surface area for: {} using mask {}".format(emmap_path.split("/")[-1], mask_path.split("/")[-1]))
    emmap = mrcfile.open(emmap_path).data
    apix_tuple = tuple(mrcfile.open(emmap_path).voxel_size.tolist())
    apix = apix_tuple[0]
    origin = mrcfile.open(emmap_path).header.origin.tolist()
    
    mask = mrcfile.open(mask_path).data
    mask = (mask==1).astype(np.int_)
    
    if mask_emmap:
        emmap = emmap * mask
    
    mask_volume = mask.sum() * apix**3
    reference_mask_volume = mask_volume * 0.2  ## Thresholded at 20% of molecular volume  
    
    print("Finding reference threshold corresponding to 20% of molecular volume determined from mask = {} ang cubed ".format(reference_mask_volume))
    reference_threshold = find_volume_matching_threshold(emmap, reference_mask_volume, apix)
    print("Reference threshold found to be {:.2f}".format(reference_threshold))
    
    surface_area_emmap_at_reference_threshold = calculate_surface_area_at_threshold(
        emmap, reference_threshold=reference_threshold, apix_tuple=apix_tuple, origin=origin)
    
    num_distinct_regions_at_reference_threshold = count_distinct_regions(emmap, reference_threshold)
    
    if num_distinct_regions_at_reference_threshold > 0:
        unit_surface_area = surface_area_emmap_at_reference_threshold / num_distinct_regions_at_reference_threshold
        print("Unit surface area for {} found to be {:.2f} nm squared per discontinous region".format(emmap_path, unit_surface_area / 100))
        print("Number of discontinous region: {}".format(num_distinct_regions_at_reference_threshold))
        return unit_surface_area
    else:
        print("Error calculating unit surface are for {}: number of distinct regions is < 1")
        return 0
    
def calculate_adjusted_surface_area(emmap_path,  fsc_resolution, mask_path, b_blurrest=-300, b_sharpest=100, mask_emmap=True):
    from locscale.include.emmer.ndimage.map_tools import sharpen_maps
    from locscale.include.emmer.ndimage.profile_tools import compute_radial_profile, frequency_array, plot_radial_profile, estimate_bfactor_through_pwlf
    from locscale.include.emmer.pdb.pdb_tools import find_wilson_cutoff
    import mrcfile
    
    print("Calculating adjusted surface area for: {} using mask {}".format(emmap_path.split("/")[-1], mask_path.split("/")[-1]))
    emmap = mrcfile.open(emmap_path).data
    apix_tuple = tuple(mrcfile.open(emmap_path).voxel_size.tolist())
    apix = apix_tuple[0]
    origin = mrcfile.open(emmap_path).header.origin.tolist()
    
    mask = mrcfile.open(mask_path).data
    mask = (mask==1).astype(np.int_)
    
    if mask_emmap:
        emmap = emmap * mask
    
    mask_volume = mask.sum() * apix**3
    reference_mask_volume = mask_volume * 0.2  ## Thresholded at 20% of molecular volume  
    
    print("Finding reference threshold corresponding to 20% of molecular volume determined from mask = {} ang cubed ".format(reference_mask_volume))
    reference_threshold = find_volume_matching_threshold(emmap, reference_mask_volume, apix)
    print("Reference threshold found to be {:.2f}".format(reference_threshold))

    wilson_cutoff = find_wilson_cutoff(mask=mask, apix=apix)
    fsc_cutoff = fsc_resolution
    
    rp_emmap = compute_radial_profile(emmap)
    freq = frequency_array(rp_emmap, apix=apix)
    
    current_bfactor = estimate_bfactor_through_pwlf(freq, rp_emmap, wilson_cutoff, fsc_cutoff)[0]
    
    sharpening_bfactor = b_sharpest - current_bfactor 
    blurring_bfactor = b_blurrest - current_bfactor
    
    sharpest_map = apply_b_sharpen(emmap, apix, b_sharpen=sharpening_bfactor, d_cut = 1/fsc_cutoff)
    most_blurred_map = apply_b_sharpen(emmap, apix, b_sharpen=blurring_bfactor, d_cut = 1/fsc_cutoff)
    
    bfactor_sharpest_map = estimate_bfactor_through_pwlf(freq, compute_radial_profile(sharpest_map), wilson_cutoff, fsc_cutoff)
    bfactor_blurred_map = estimate_bfactor_through_pwlf(freq, compute_radial_profile(most_blurred_map), wilson_cutoff, fsc_cutoff)
    
    print("Current bfactor \t: {}".format(current_bfactor))
    print("Bfactor of sharpest map \t: {}".format(bfactor_sharpest_map[0]))
    print("Bfactor of blurred map \t: {}".format(bfactor_blurred_map[0]))
    
    scale_factor, asa_blurred, asa_sharpened = find_c_scale(
        sharpest_map, most_blurred_map, reference_threshold=reference_threshold, apix_tuple=apix_tuple, origin=origin)
    
    surface_area_emmap_at_reference_threshold = calculate_surface_area_at_threshold(
        emmap, reference_threshold=reference_threshold, apix_tuple=apix_tuple, origin=origin)
    
    num_distinct_regions_at_reference_threshold = count_distinct_regions(emmap, reference_threshold)
    
    adjusted_surface_area = surface_area_emmap_at_reference_threshold - scale_factor * num_distinct_regions_at_reference_threshold
    
    print("Adjusted surface area measured to be: {:.2f} nm squared".format(adjusted_surface_area/100))
    
    if adjusted_surface_area < asa_blurred or adjusted_surface_area < asa_sharpened:
        print("Calculated adjusted surface area lesser than the extreme values of bfactor considered for calculated scale factor. ")
        print("Adjusted surface area at the most blurred: {:.2f} nm squared".format(asa_blurred/100))
        print("Adjusted surface area at the most sharpened: {:.2f} nm squared".format(asa_sharpened/100))
    
    return adjusted_surface_area


#%%


folder = "/mnt/c/Users/abharadwaj1/Downloads/ForUbuntu/LocScale/tests/map_quality/emd5778"

emmap_path = os.path.join(folder, "loc_scale_oct14.mrc")
mask_path = os.path.join(folder, "pdb3j5p_mask.mrc")
fsc_resolution = 3.4
    

#%%
