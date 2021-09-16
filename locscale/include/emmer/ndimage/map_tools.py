#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Program to symmetrize maps
Created on Thu Apr 15 00:40:10 2021

@author: alok
"""
import numpy as np
import mrcfile
import gemmi


from locscale.include.emmer.ndimage.map_utils import parse_input
    
def compute_real_space_correlation(input_map_1,input_map_2):
    '''
    Function to calculate the Real Space Cross Correlation (RSCC) between two maps, or any two ndarrays. 
    
    RSCC is calculated by standardizing two arrays by subtracting their mean and dividing by their standard deviation

    Parameters
    ----------
    array1 : numpy.ndarray
        
    array2 : numpy.ndarray
        

    Returns
    -------
    RSCC : float
        Floating point number between 0 and 1 showing the RSCC between two arrays

    '''
    from locscale.include.emmer.ndimage.map_utils import parse_input
    array1 = parse_input(input_map_1, allow_any_dims=True)
    array2 = parse_input(input_map_2, allow_any_dims=True)
    
    (map1_mean,map1_std) = (array1.mean(),array1.std())
    (map2_mean,map2_std) = (array2.mean(),array2.std())
    
    n = array1.size
    
    RSCC = (((array1-map1_mean)*(array2-map2_mean))/(map1_std*map2_std)).sum() * (1/n)
    
    return RSCC

def get_center_of_mass(emmap_data, apix):
    '''
    Computes the center of mass of a given input emmap. 
    Note: converts the negative intensities to positive to calculate COM

    Parameters
    ----------
    emmap_data : numpy.ndarray
        
    apix : float or any iterable
        Voxelsize

    Returns
    -------
    com_real : numpy.ndarray
        units: (A * A * A)

    '''
    from scipy.ndimage import center_of_mass
    from locscale.include.emmer.ndimage.map_utils import convert_to_tuple
    
    com_pixels = np.array(center_of_mass(abs(emmap_data)))
    apix_array = np.array(convert_to_tuple(apix))
    
    com_real = com_pixels * apix_array
    
    return com_real
    
def add_half_maps(halfmap_1_path, halfmap_2_path, output_filename):
    '''
    Function to add two half maps

    Parameters
    ----------
    halfmap_1_path : TYPE
        DESCRIPTION.
    halfmap_2_path : TYPE
        DESCRIPTION.

    Returns
    -------
    None.

    '''
    import mrcfile
    from locscale.include.emmer.ndimage.map_utils import save_as_mrc
    halfmap1 = mrcfile.open(halfmap_1_path).data
    halfmap2 = mrcfile.open(halfmap_2_path).data
    
    full_map = halfmap1 + halfmap2
    full_voxel_size = mrcfile.open(halfmap_1_path).voxel_size
    full_header = mrcfile.open(halfmap_1_path).header
    save_as_mrc(map_data=full_map, output_filename=output_filename, apix=full_voxel_size, verbose=True, header=full_header) 
    
    return output_filename

def add_half_maps_2(halfmap_1_path, halfmap_2_path, output_filename):
    '''
    Function to add two half maps

    Parameters
    ----------
    halfmap_1_path : TYPE
        DESCRIPTION.
    halfmap_2_path : TYPE
        DESCRIPTION.

    Returns
    -------
    None.

    '''
    import mrcfile
    from locscale.include.emmer.ndimage.map_utils import save_as_mrc_2
    halfmap1 = mrcfile.open(halfmap_1_path).data
    halfmap2 = mrcfile.open(halfmap_2_path).data
    
    full_map = halfmap1 + halfmap2
    full_voxel_size = mrcfile.open(halfmap_1_path).voxel_size
    full_header = mrcfile.open(halfmap_1_path).header
    save_as_mrc_2(map_data=full_map, output_filename=output_filename, apix=full_voxel_size, verbose=True, header=full_header) 
    
    return output_filename
    
def estimate_global_bfactor_map(emmap, apix, wilson_cutoff, fsc_cutoff):
    from locscale.include.emmer.ndimage.profile_tools import number_of_segments, frequency_array, compute_radial_profile, estimate_bfactor_through_pwlf
    
    rp_unsharp = compute_radial_profile(emmap)
    freq = frequency_array(amplitudes=rp_unsharp, apix=apix)
    num_segments = number_of_segments(fsc_cutoff)
        
    bfactor,_,(fit,z,slopes) = estimate_bfactor_through_pwlf(freq,rp_unsharp, wilson_cutoff=wilson_cutoff, fsc_cutoff=fsc_cutoff, num_segments=num_segments)
    
    return bfactor, z, slopes, fit
    
def compute_scale_factors(em_profile, ref_profile):
    scale_factor = np.sqrt(ref_profile**2/em_profile**2)
    return scale_factor

def set_radial_profile(vol, scale_factor, radii):
    ps = np.fft.rfftn(vol)
    for j,r in enumerate(np.unique(radii)[0:vol.shape[0]//2]):
            idx = radii == r
            ps[idx] *= scale_factor[j]

    return np.fft.irfftn(ps, s=vol.shape)    

def sharpen_maps(vol, apix, global_bfactor=0):
    '''
    Function to apply a global sharpening factor to EM density maps 

    Parameters
    ----------
    vol : numpy.ndarray (dims=3)
        Input map
    apix : Float
        Pixelsize (one dimension only)
    global_bfactor : int, optional
        The default is 0.

    Returns
    -------
    sharpened_map : numpy.ndarray (dims=3)

    '''
    from locscale.include.emmer.ndimage.profile_tools import compute_radial_profile, frequency_array
    from locscale.include.emmer.ndimage.map_tools import set_radial_profile, compute_scale_factors
    
    emmap_profile, radii = compute_radial_profile(vol,return_indices=True)
    freq = frequency_array(amplitudes=emmap_profile, apix=apix)
    
    sharpened_profile = emmap_profile * np.exp(-1*global_bfactor/4 * freq**2)

    scale_factors = compute_scale_factors(emmap_profile, sharpened_profile)
    sharpened_map = set_radial_profile(vol, scale_factors, radii)
    
    return sharpened_map
    
def crop_map_between_residues(emmap_path, pdb_path, chain_name, residue_range=None, dilation_radius=3):
    '''
    Function to extract map intensities around atoms between a given residue range

    Parameters
    ----------
    emmap_path : str
        Path to a map file 
    pdb_path : str
        Path to a PDB/MMCIF file
    chain_name : str
        Chain name
    residue_range : list, optional
        To extract all atoms between residue id
        residue_range=[start_res_id, end_res_id] (both incl). The default is [0,-1], which returns all residues present in a chain. 
    dilation_radius : float, optional
        The radius of the sphere (in Ang) to place at atomic positions determined by the PDB file. Default is 3A.

    Returns
    -------
    cropped_map : numpy.ndarray
    

    '''
    from locscale.include.emmer.pdb.pdb_tools import get_atomic_positions_between_residues
    from locscale.include.emmer.ndimage.map_utils import convert_pdb_to_mrc_position, dilate_mask
    
    apix = mrcfile.open(emmap_path).voxel_size.x
    emmap = mrcfile.open(emmap_path).data
    
    map_shape = emmap.shape
    
    mask = np.zeros(map_shape)
    
    gemmi_st = gemmi.read_structure(pdb_path)
    
    pdb_positions = get_atomic_positions_between_residues(gemmi_st, chain_name, residue_range)
    
    print("Found {} atom sites".format(len(pdb_positions)))
    
    mrc_position = convert_pdb_to_mrc_position(pdb_positions, apix)
    zz,yy,xx = zip(*mrc_position)
    mask[zz,yy,xx] = 1
    
    #dilation_radius = 3 #A
    dilation_radius_int = round(dilation_radius / apix)
    dilated_mask = dilate_mask(mask, radius=dilation_radius_int)
    
    cropped_map = emmap * dilated_mask
    
    return cropped_map

    