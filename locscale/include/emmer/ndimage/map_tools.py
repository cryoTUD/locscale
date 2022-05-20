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
    halfmap_1_path : str
        
    halfmap_2_path : str
        

    Returns
    -------
    output_filename : str

    '''
    import mrcfile
    from locscale.include.emmer.ndimage.map_utils import save_as_mrc
    halfmap1 = mrcfile.open(halfmap_1_path).data
    halfmap2 = mrcfile.open(halfmap_2_path).data
    
    if halfmap1.shape == halfmap2.shape:
       full_map = (halfmap1 + halfmap2)/2
       full_voxel_size = mrcfile.open(halfmap_1_path).voxel_size.tolist()
    
       full_header = mrcfile.open(halfmap_1_path).header
       #print(full_voxel_size)
       save_as_mrc(map_data=full_map, output_filename=output_filename, apix=full_voxel_size, verbose=True) 
    
       return output_filename
    else:
      print("Half maps are not of equal dimension.")

    
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


def compute_radial_profile_proper(vol, frequency_map):

    vol_fft = np.fft.rfftn(vol, norm="ortho");
    dim = vol_fft.shape;
    ps = np.real(np.abs(vol_fft));
    frequencies = np.fft.rfftfreq(dim[0]);
    #bins = np.digitize(frequency_map, frequencies);
    #bins = bins - 1;
    x, y, z = np.indices(ps.shape)
    radii = np.sqrt(x**2 + y**2 + z**2)
    radii = radii.astype(int)
    radial_profile = np.bincount(radii.ravel(), ps.ravel()) / np.bincount(radii.ravel())
    radial_profile = radial_profile[0:int(ps.shape[0]/2)+1]

    return radial_profile, frequencies;

def compute_radial_profile_simple(vol, return_frequencies=False):
    from locscale.include.confidenceMapUtil import FDRutil
    frequency_map = FDRutil.calculate_frequency_map(np.zeros(vol.shape))
    
    em_profile, frequencies_map = compute_radial_profile_proper(vol, frequency_map)
    
    if return_frequencies:
        return em_profile, frequencies_map
    else:
        return em_profile

def set_radial_profile_simple(vol, scale_factors, frequencies):
    from locscale.include.confidenceMapUtil import FDRutil
    frequency_map = FDRutil.calculate_frequency_map(np.zeros(vol.shape))
    
    map_shape = vol.shape
    map_b_sharpened, _ = set_radial_profile_proper(vol, scale_factors, frequencies, frequency_map, map_shape);
    
    return map_b_sharpened
    
    

def set_radial_profile_proper(vol, scale_factors, frequencies, frequency_map, shape):
    vol_fft = np.fft.rfftn(np.copy(vol), norm='ortho');
    scaling_map = np.interp(frequency_map, frequencies, scale_factors);
    scaled_map_fft = scaling_map * vol_fft;
    scaled_map = np.real(np.fft.irfftn(scaled_map_fft, shape, norm='ortho'));

    return scaled_map, scaled_map_fft;

def set_radial_profile(vol, scale_factor, radii):
    ps = np.fft.rfftn(vol)
    for j,r in enumerate(np.unique(radii)[0:vol.shape[0]//2]):
            idx = radii == r
            ps[idx] *= scale_factor[j]

    return np.fft.irfftn(ps, s=vol.shape)  

def set_radial_profile_to_volume(emmap, ref_profile):
    from locscale.include.emmer.ndimage.profile_tools import compute_radial_profile
    em_profile,frequencies = compute_radial_profile_simple(emmap, return_frequencies=True)
    scale_factors = compute_scale_factors(em_profile, ref_profile)
    scaled_map = set_radial_profile_simple(emmap, scale_factors, frequencies)
    
    return scaled_map

    

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

def get_atomic_model_mask(emmap_path, pdb_path, dilation_radius=3, softening_parameter=5, output_filename=None):
    '''
    Function to extract map intensities around atoms between a given residue range

    Parameters
    ----------
    emmap_path : str
        Path to a reference emmap for metadata
    pdb_path : str
        Path to a PDB/MMCIF file
  
    dilation_radius : float, optional
        The radius of the sphere (in Ang) to place at atomic positions determined by the PDB file. Default is 3A.

    Returns
    -------
    model_mask : str
    

    '''
    from locscale.include.emmer.pdb.pdb_tools import get_atomic_positions_between_residues
    from locscale.include.emmer.ndimage.map_utils import convert_pdb_to_mrc_position, dilate_mask, save_as_mrc
    import mrcfile
    import os
    from locscale.include.emmer.ndimage.filter import get_cosine_mask
    
    pdb_folder = os.path.dirname(pdb_path)
    pdb_name = os.path.basename(pdb_path)
    apix = mrcfile.open(emmap_path).voxel_size.tolist()[0]
    emmap = mrcfile.open(emmap_path).data
    map_shape = emmap.shape
    
    
    gemmi_st = gemmi.read_structure(pdb_path)
    
    
    mask = np.zeros(map_shape)
    pdb_positions = []
    for model in gemmi_st:
        for chain in model:
            for res in chain:
                for atom in res:
                    pdb_positions.append(atom.pos.tolist())
                        
        
    print("Found {} atom sites".format(len(pdb_positions)))
        
    mrc_position = convert_pdb_to_mrc_position(pdb_positions, apix)
    zz,yy,xx = zip(*mrc_position)
    mask[zz,yy,xx] = 1
        
    dilation_radius_int = round(dilation_radius / apix)
    dilated_mask = dilate_mask(mask, radius=dilation_radius_int)
    
    softened_mask = get_cosine_mask(dilated_mask, length_cosine_mask_1d=softening_parameter)
    

    if output_filename is None:
        output_filename = os.path.join(pdb_folder, pdb_name[:-4]+"_model_mask.mrc")
            
    save_as_mrc(softened_mask, output_filename=output_filename, apix=apix)
        
    return output_filename
    
    
def apply_radial_profile(emmap, reference_map):
    from locscale.include.emmer.ndimage.map_tools import set_radial_profile, compute_scale_factors
    from locscale.include.emmer.ndimage.profile_tools import compute_radial_profile, frequency_array
    from locscale.include.emmer.ndimage.map_utils import parse_input
    
    emmap = parse_input(emmap)
    reference_map = parse_input(reference_map)
    
    rp_emmap = compute_radial_profile(emmap)
    rp_reference_map, radii = compute_radial_profile(reference_map, return_indices=True)
    
    sf = compute_scale_factors(rp_emmap, rp_reference_map)
    
    scaled_map = set_radial_profile(emmap, sf, radii)
    
    return scaled_map

def get_local_bfactor_emmap(emmap_path, center, fsc_resolution, boxsize=None, standard_notation=True, mask_path=None, wilson_cutoff="singer"):
    from locscale.include.emmer.ndimage.profile_tools import estimate_bfactor_standard, compute_radial_profile, frequency_array, plot_radial_profile
    from locscale.include.emmer.ndimage.map_tools import compute_real_space_correlation
    from locscale.include.emmer.ndimage.map_utils import measure_mask_parameters, get_all_voxels_inside_mask, extract_window
    from locscale.include.emmer.pdb.pdb_tools import find_wilson_cutoff
    from locscale.utils.math_tools import round_up_to_even
    import random
    import mrcfile
    from tqdm import tqdm

    emmap = mrcfile.open(emmap_path).data
    apix = mrcfile.open(emmap_path).voxel_size.tolist()[0]
    if wilson_cutoff == "singer":
        mask = mrcfile.open(mask_path).data
        global_wilson_cutoff = find_wilson_cutoff(mask_path=mask_path, verbose=False)
    else:
        global_wilson_cutoff = 10
    
    fsc_cutoff = fsc_resolution
    if boxsize is None:
        boxsize = round_up_to_even(25 / apix)
    else:
        boxsize = round_up_to_even(boxsize)

    
    
    
            
    emmap_window = extract_window(emmap, center, boxsize)
            
    rp_local = compute_radial_profile(emmap_window)
    freq = frequency_array(rp_local, apix)
    if wilson_cutoff == "singer":                       
        mask_window = extract_window(mask, center, boxsize)
        num_atoms, _  = measure_mask_parameters(mask=mask_window, apix=apix, verbose=False)
        local_wilson_cutoff = find_wilson_cutoff(num_atoms=num_atoms)
        local_wilson_cutoff = np.clip(local_wilson_cutoff, fsc_cutoff*1.5, global_wilson_cutoff)
    else:
        local_wilson_cutoff = 10
            
            
    bfactor,qfit = estimate_bfactor_standard(freq, rp_local, local_wilson_cutoff, fsc_cutoff, standard_notation=standard_notation, return_fit_quality=True)
            
    return bfactor, qfit

def get_bfactor_distribution(emmap_path, mask_path, fsc_resolution, boxsize=None, num_centers=15000, standard_notation=True, wilson_cutoff="singer"):
    from locscale.include.emmer.ndimage.profile_tools import estimate_bfactor_standard, compute_radial_profile, frequency_array, plot_radial_profile
    from locscale.include.emmer.ndimage.map_tools import compute_real_space_correlation
    from locscale.include.emmer.ndimage.map_utils import measure_mask_parameters, get_all_voxels_inside_mask, extract_window
    from locscale.include.emmer.pdb.pdb_tools import find_wilson_cutoff
    from locscale.utils.math_tools import round_up_to_even
    import random
    import mrcfile
    from tqdm import tqdm

    emmap = mrcfile.open(emmap_path).data
    mask = mrcfile.open(mask_path).data
    
    apix = mrcfile.open(mask_path).voxel_size.tolist()[0]
    global_wilson_cutoff = find_wilson_cutoff(mask_path=mask_path)
    fsc_cutoff = fsc_resolution
    if boxsize is None:
        boxsize = round_up_to_even(25 / apix)
    else:
        boxsize = round_up_to_even(boxsize)
    print(boxsize)
    all_points = get_all_voxels_inside_mask(mask_input=mask, mask_threshold=1)
    random_centers = random.sample(all_points,num_centers)
    
    bfactor_distributions = {}
    
    for center in tqdm(random_centers, desc="Analysing local bfactors distribution"):
        try:
            
            emmap_window = extract_window(emmap, center, boxsize)
            
            rp_local = compute_radial_profile(emmap_window)
            freq = frequency_array(rp_local, apix)
            if wilson_cutoff == "singer":                       
                mask_window = extract_window(mask, center, boxsize)
                num_atoms, _  = measure_mask_parameters(mask=mask_window, apix=apix, verbose=False)
                local_wilson_cutoff = find_wilson_cutoff(num_atoms=num_atoms)
                local_wilson_cutoff = np.clip(local_wilson_cutoff, fsc_cutoff*1.5, global_wilson_cutoff)
            else:
                local_wilson_cutoff = 10
            
            
            bfactor,qfit = estimate_bfactor_standard(freq, rp_local, local_wilson_cutoff, fsc_cutoff, standard_notation=standard_notation, return_fit_quality=True)
            
            bfactor_distributions[tuple(center)] = tuple([bfactor, qfit])
        except Exception as e:
            print("Error at {}".format(center))
            print(e)
            raise
            
        
    
    return bfactor_distributions

def get_bfactor_distribution_multiple(list_of_emmap_paths, mask_path, fsc_resolution, boxsize=None, num_centers=15000, standard_notation=True, wilson_cutoff="singer"):
    from locscale.include.emmer.ndimage.profile_tools import estimate_bfactor_standard, compute_radial_profile, frequency_array, plot_radial_profile
    from locscale.include.emmer.ndimage.map_tools import compute_real_space_correlation
    from locscale.include.emmer.ndimage.map_utils import measure_mask_parameters, get_all_voxels_inside_mask, extract_window
    from locscale.include.emmer.pdb.pdb_tools import find_wilson_cutoff
    from locscale.utils.math_tools import round_up_to_even
    import random
    import mrcfile
    from tqdm import tqdm
    import os

    
    mask = mrcfile.open(mask_path).data
    
    apix = mrcfile.open(mask_path).voxel_size.tolist()[0]
    global_wilson_cutoff = find_wilson_cutoff(mask_path=mask_path)
    fsc_cutoff = fsc_resolution
    if boxsize is None:
        boxsize = round_up_to_even(25 / apix)
    else:
        boxsize = round_up_to_even(boxsize)
    print(boxsize)
    all_points = get_all_voxels_inside_mask(mask_input=mask, mask_threshold=1)
    random_centers = random.sample(all_points,num_centers)
    
    bfactor_distributions = {}
    
    for emmap_path in list_of_emmap_paths:
        emmap_name = os.path.basename(emmap_path)
        emmap = mrcfile.open(emmap_path).data
        temp_distribution = {}
        for center in tqdm(random_centers, desc="Analysing local bfactors distribution"):
            
            try:
                
                emmap_window = extract_window(emmap, center, boxsize)
                
                rp_local = compute_radial_profile(emmap_window)
                freq = frequency_array(rp_local, apix)
                if wilson_cutoff == "singer":                       
                    mask_window = extract_window(mask, center, boxsize)
                    num_atoms, _  = measure_mask_parameters(mask=mask_window, apix=apix, verbose=False)
                    local_wilson_cutoff = find_wilson_cutoff(num_atoms=num_atoms)
                    local_wilson_cutoff = np.clip(local_wilson_cutoff, fsc_cutoff*1.5, global_wilson_cutoff)
                else:
                    local_wilson_cutoff = 10
                
                
                bfactor,qfit = estimate_bfactor_standard(freq, rp_local, local_wilson_cutoff, fsc_cutoff, standard_notation=standard_notation, return_fit_quality=True)
                
                temp_distribution[tuple(center)] = tuple([bfactor, qfit])
            except Exception as e:
                print("Error at {}".format(center))
                print(e)
                raise
        bfactor_distributions[emmap_name] = temp_distribution
            
        
    
    return bfactor_distributions
    
        
        
