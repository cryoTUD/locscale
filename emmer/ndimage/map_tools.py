#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Program to symmetrize maps
Created on Thu Apr 15 00:40:10 2021

@author: alok
"""
from emmer.headers import  *
from emmer.ndimage.map_utils import parse_input
    
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
    from emmer.ndimage.map_utils import parse_input
    array1 = parse_input(input_map_1)
    array2 = parse_input(input_map_2)
    
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
    from emmer.ndimage.map_utils import convert_to_tuple
    
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
    from emmer.ndimage.map_utils import save_as_mrc
    halfmap1 = mrcfile.open(halfmap_1_path).data
    halfmap2 = mrcfile.open(halfmap_2_path).data
    
    full_map = halfmap1 + halfmap2
    voxel_size = mrcfile.open(halfmap_1_path).voxel_size
    save_as_mrc(map_data=full_map, output_filename=output_filename, header=mrcfile.open(halfmap_1_path).header)
    
    return output_filename
    
def symmetrize(emmap_path, res, pg, save_mrc=True, voxel_size=None):
    '''
    Function to return a symmtery averaged map given an input map. 

    Parameters
    ----------
    emmap_path : string
        Can input either a path to MRC file or the 
    
    res : float
        Resolution of the structure in Angstrom
    pg : string
        Point groups symmetry ('C3, 'O', etc)
    save_mrc : TYPE, optional
        DESCRIPTION. The default is True.
        
    

    Returns
    -------
    If save_mrc flag is set true, then the function saves the symmetry averaged map
    and returns the filepath and ndarray. Else, it only returns the ndarray

    '''
    import emda_methods as em #pythonpath contains  path/to/emda/emda/ 
    
    sym = em.symmetry_average([emmap_path],[res],pglist=[pg])
    apix = mrcfile.open(emmap_path).voxel_size
	
    if os.path.exists(filename) and save_mrc:
        print("Applied symmetry averaging successfully. Find the map here: \n"+filename+"\nRSCC with input map = "+str(compute_real_space_correlation(mrcfile.open(emmap_path).data, mrcfile.open(filename).data)))
        filename=emmap_path[:-4]+"_"+pg+"_symmetry.mrc"
        save_as_mrc(sym[0], apix, filename)
        return filename, sym[0]
    elif os.path.exists(filename) and not save_mrc:
        return sym[0]
    

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
    from emmer.ndimage.profile_tools import compute_radial_profile, frequency_array
    from emmer.ndimage.map_tools import set_radial_profile, compute_scale_factors
    
    emmap_profile, radii = compute_radial_profile(vol,return_indices=True)
    freq = frequency_array(amplitudes=emmap_profile, apix=apix)
    
    sharpened_profile = emmap_profile * np.exp(-1*global_bfactor/4 * freq**2)

    scale_factors = compute_scale_factors(emmap_profile, sharpened_profile)
    sharpened_map = set_radial_profile(vol, scale_factors, radii)
    
    return sharpened_map
    
def crop_map_between_residues(emmap_path, pdb_path, chain_name, residue_range=[0,-1], dilation_radius=3):
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
    from emmer.pdb.pdb_tools import get_atomic_positions_between_residues
    from emmer.ndimage.util import convert_pdb_to_mrc_position, dilate_mask
    
    apix = mrcfile.open(emmap_path).voxel_size.x
    emmap = mrcfile.open(emmap_path).data
    
    map_shape = emmap.shape
    
    mask = np.zeros(map_shape)
    
    gemmi_st = gemmi.read_structure(pdb_path)
    
    pdb_positions = get_atomic_positions_between_residues(gemmi_st, chain_name, residue_range)
    
    mrc_position = convert_pdb_to_mrc_position(pdb_positions, apix)
    zz,yy,xx = zip(*mrc_position)
    mask[zz,yy,xx] = 1
    
    #dilation_radius = 3 #A
    dilation_radius_int = round(dilation_radius / apix)
    dilated_mask = dilate_mask(mask, radius=dilation_radius_int)
    
    cropped_map = emmap * dilated_mask
    
    return cropped_map

    