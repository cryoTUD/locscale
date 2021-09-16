#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec  3 14:11:29 2020

@author: Alok, Arjen, Maarten, Stefan, Reinier 
"""

# pdb_tools contains several useful fucntions for common manipulations with
# pdb structures, making use of the gemmi package. functions are classified 
# into pdb_tools when they can be considered an application on their own
# but do not have so many distinct features that they warrent their own script.

# global imports
import mrcfile
import gemmi
import numpy as np
import json
import pypdb
import os
import sys
from scipy import signal
from emmer.ndimage.filter import *
from scipy.optimize import curve_fit
from scipy.interpolate import interp1d

#%% functions

def add_cryst1_line(pdb_path,unitcell=None,emmap_path=None,new_pdb_path=None):
    '''
    pdb_path -> Address of .pdb path
    
    Some PDB files developed for cryoEM maps do not have proper cryst1 record. Two options to modify:

    1. From an input tuple, or array. In this case, unitcell is a python tuple, which has unit cell dimensions in angstorm
    Ex: unitcell = (x,y,z)
    2. From a mrcfile. In this case, point to an associated EM map and the unit cell dimensions are taken from that
    emmap_path -> Address of associated .mrc file
    
    If you like to the pdb file with a different name, or address then change the 'new_pdb_path' 
    
    '''
    if emmap_path is not None:
        mrc = mrcfile.open(emmap_path)
        cella = mrc.header.cella
        x = cella.x
        y = cella.y
        z = cella.z
    elif unitcell is not None:
        x = unitcell[0]
        y = unitcell[1]
        z = unitcell[2]
    else:
        print("Please give either unit cell dimensions (in Ang) or point to an associated mrc file!")
        return
    
    unitcell = gemmi.UnitCell(x,y,z,90,90,90)
    
    gemmi_structure = gemmi.read_structure(pdb_path)
    gemmi_structure.cell = unitcell
    if new_pdb_path is None:
        gemmi_structure.write_pdb(pdb_path)
    else:
        gemmi_structure.write_pdb(new_pdb_path)
        
def set_to_center_of_unit_cell(pdb_structure, unitcell):
    '''
    Function to set the center of mass of a PDB structure to the center of a unitcell

    Parameters
    ----------
    pdb_structure : gemmi.Structure
        Input structure 
    unitcell : gemmi.UnitCell
        Input unitcell

    Returns
    -------
    centered_pdb : gemmi.Structure

    '''
    from emmer.pdb.pdb_utils import shift_coordinates
    
    pdb_structure_local = pdb_structure.clone()
    center_of_mass_old = np.array(pdb_structure_local[0].calculate_center_of_mass().tolist())
    center_of_mass_new = np.array([unitcell.a/2, unitcell.b/2, unitcell.c/2])
    
    translation_matrix = center_of_mass_new - center_of_mass_old
    shifted_structure = shift_coordinates(trans_matrix=translation_matrix, input_structure=pdb_structure_local)
    
    return shifted_structure
    
    

def get_unit_cell_estimate(pdb_struct,vsize):
          
    '''
    Find an estimated size of unit cell in A based on nunber of atoms and apix

    As reference: PDB3J5P had ~18000 atoms and a box size of 256^3
          
    '''

    number_of_atoms = pdb_struct[0].count_atom_sites()
    estimated_box_size = number_of_atoms * 256 / 18000
    unit_cell_dim =  estimated_box_size * vsize
    unitcell = gemmi.UnitCell(unit_cell_dim,unit_cell_dim,unit_cell_dim,90,90,90)
          
    return unitcell
        
# def pdb_to_map(
#         pdb_id=None,pdb_path=None,pdb_structure=None,mdlidx=0,resolution=None,vsize=None,unitcell=None,save_mrc_path=None,
#         set_zero_origin=False,crop_to_center=False, perform_fftshift=False, verbose=True,remove_waters=True):
#      '''

#     Parameters
#     ----------
#     pdb_id : string, optional
#         PDB ID of the coordinate file, incase the file is not present locally or a gemmi structure is absent. Example: pdb_id = '3j5p'
#     pdb_path : string, optional
#         Location of .pdb file, incase the file is present locally. 
#     pdb_structure : gemmi.Structure(), Note: Either pdb_id, pdb_path or pdb_structure is required.
#         Gemmi structure file for programs written using gemmi API. 
#     resolution : float, optional
#         Resolution cut-off in Angstrom. 
#     vsize : float, required
#         apix in Angstrom. 
#     unitcell : gemmi.UnitCell(), optional
#         If this parameter is passed, then maps will be simulated strictly using this. Example: unitcell = gemmi.UnitCell(311.1,311.1,311.1,90,90,90)
#     save_mrc_path : string, optional
#         If this parameter is passed, then maps simulated will be stored locally at this location. 
#         Please provide the .mrc extension. Example: "/dev/test.mrc"
#     set_zero_origin : BOOL, optional
#         If this option is set True, then the coordinates will be shifted such that center of mass will be at zero origin. 
#     perform_fftshift : BOOL, optional
#         If this option is set True, then  np.fft.fftshift(,axes=(0,1,2)) will be performed for the simulated volume
#     verbose : BOOL, optional
#         If set to false, the print statements in the function are suppressed.

#     Returns
#     -------
#     emmap : numpy.ndarray
#         This is the simulated map showing the 3D volume of the structure

#     '''

#      if pdb_id is not None:
#           pdbfile = pypdb.get_pdb_file(pdb_id,filetype='pdb',compression=False)
#           pdb_struct = gemmi.read_pdb_string(pdbfile)
          
#           # model = pdb_struct[0]
#      elif pdb_path is not None:
#           pdb_struct = gemmi.read_pdb(pdb_path)
          
#           # model = pdb_struct[0]
#      elif pdb_structure is not None:
#           pdb_struct = pdb_structure
          
#           # model = pdb_struct[0]
#      else:
#           print("Please input a PDB model or path to PDB file or a PDB ID (as string)!")
#           return 0
     
#      if abs(pdb_struct.cell.a * pdb_struct.cell.b * pdb_struct.cell.c) <= 64 and unitcell is None:
#           pdb_struct.cell = get_unit_cell_estimate(pdb_struct,vsize)
          
#      elif unitcell is not None:
#           pdb_struct.cell = unitcell
          
#      if vsize is None:
#          print("Please enter a voxelsize (in A). Returning null")
#          return 0

        
#      if mdlidx < len(pdb_struct):
#         model = pdb_struct[mdlidx]
#      model_origin = np.array(model.calculate_center_of_mass().tolist()) # Center of mass of protein is taken as the "model_origin". This should be set to zero
     
#      # Check if set_zero_origin is True and if model_origin is more than 1A away from zero origin
#      if set_zero_origin:
         
#          expected_side_length_ang = pdb_struct.cell.a
#          expected_side_length_pix = expected_side_length_ang / vsize
         
#          new_origin = np.array([expected_side_length_pix//2,expected_side_length_pix//2,expected_side_length_pix//2])
         
#          if verbose:
#               print("Setting model coordinates to zero origin")
#               print("Initial Origin: "+str(model_origin.round(2)))
#          translation_matrix = new_origin - model_origin
#          model = shift_coordinates(trans_matrix=translation_matrix, input_model=model)
#          new_model_origin = np.array(model.calculate_center_of_mass().tolist()) 
#          if verbose:
#              print("Final Origin: "+str(new_model_origin.round(2)))
     
#      if verbose:
#          print("Using Gemmi structure with unit cell dimensions",pdb_struct.cell)
     
#      if remove_waters:
#          pdb_struct.remove_waters()
#      pdb_struct.remove_empty_chains()
     
#      dencalc = gemmi.DensityCalculatorE()
#      dencalc.d_min = vsize * 2  # Added the weird coefficient at the end to get expected size in pixel (Can change, need to check!)
#      dencalc.rate = 1
#      dencalc.blur= 0
#      dencalc.set_grid_cell_and_spacegroup(pdb_struct)
#      dencalc.put_model_density_on_grid(model)
       
#      emmap = np.array(dencalc.grid,copy=False)
     
#      if perform_fftshift:
#          emmap = np.fft.fftshift(emmap)
         
#      if resolution is not None:
#          emmap = low_pass_filter(emmap,resolution,vsize)
          
#      if crop_to_center:
#          X,Y,Z = emmap.shape
#          cx,cy,cz = X//2,Y//2,Z//2
         
#          expected_half_side_length = int(expected_side_length_pix//2)
#          emmap = emmap[cz-expected_half_side_length:cz+expected_half_side_length,cy-expected_half_side_length:cy+expected_half_side_length,cx-expected_half_side_length:cx+expected_half_side_length]
#          emmap = np.rot90(emmap,k=3,axes=(0,2))
          
#      if save_mrc_path is not None:
#           voxel_size_record = np.rec.array((vsize,vsize,vsize),dtype=[('x','<f4'),('y','<f4'),('z','<f4')])
#           save_as_mrc(emmap,voxel_size_record,save_mrc_path)       
          
#      return emmap

def find_radius_of_gyration(model_path=None, input_gemmi_st=None):
    if model_path is not None:
        gemmi_st = gemmi.read_pdb(model_path)
    elif input_gemmi_st is not None:
        gemmi_st = input_gemmi_st.clone()
    else:
        print("Input error!")
        return 0
    
    num_atoms = gemmi_st[0].count_atom_sites()
    com = gemmi_st[0].calculate_center_of_mass()
    distances = []
    for model in gemmi_st:
        for chain in model:
            for res in chain:
                ca = res.get_ca()
                if ca is not None:
                    distances.append(com.dist(ca.pos))
    
    np_distance = np.array(distances)
    
    Rg = np.sum(np_distance**2)/num_atoms
    
    return Rg

def find_wilson_cutoff(model_path=None, input_gemmi_st=None, mask_path=None, mask = None, apix=None, method='Singer', return_as_frequency=False):
    '''
    Function to find the cutoff frequency above which Wilson statistics hold true. If a PDB file is passed either as a gemmi structure as a PDB path, then radius of gyration is found rigorously by the mean distance to center of mass of protein. If a mask is passed, however, then radius of gyration is estimated from the num_atoms calculated from the mask volume. 
    
Reference: 
    1) Estimating Radius of gyration from num_atoms John J. Tanner,  2016 (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5053138/)
    
    2) Estimating cutoff frequency: Amit Singer, 2021 (https://www.biorxiv.org/content/10.1101/2021.05.14.444177v1.full)
    
    3) Estimating cutoff frequency: Guiner method - Rosenthal & Henderson, 2003 (https://doi.org/10.1016/j.jmb.2003.07.013)

    Parameters
    ----------
    model_path : string, optional
        path to pdb file. The default is None.
    input_gemmi_st : gemmi.Structure(), optional
        
    mask_path : string, optional
        path to mask. The default is None.
    method : string, optional
        Method used to find the cutoff frequency. Two accepted values are: 'Singer', and 'Rosenthal_Henderson' (case insensitive). The default is 'Singer'.
    return_as_frequency : TYPE, optional
        DESCRIPTION. The default is False.

    Returns
    -------
    TYPE
        DESCRIPTION.

    '''
    
    from emmer.ndimage.map_utils import measure_mask_parameters
    if model_path is not None:
        gemmi_st = gemmi.read_pdb(model_path)
        num_atoms = gemmi_st[0].count_atom_sites()
        Rg = find_radius_of_gyration(input_gemmi_st=gemmi_st)
    elif input_gemmi_st is not None:
        gemmi_st = input_gemmi_st.clone()
        num_atoms = gemmi_st[0].count_atom_sites()
        Rg = find_radius_of_gyration(input_gemmi_st=gemmi_st)
    elif mask_path is not None:
        mask_vol_A3, protein_mass, num_atoms, mask_dims,maskshape = measure_mask_parameters(mask_path=mask_path, detailed_report=True)
        
        R_constant = 2 #A
        v = 0.4 # Exponent derived empirically Ref. 1 for monomers and oligomers
        Rg = R_constant * num_atoms**v
    elif mask is not None and apix is not None and mask_path is None:
        mask_vol_A3, protein_mass, num_atoms, mask_dims,maskshape = measure_mask_parameters(mask=mask, apix=apix, detailed_report=True)
        
        R_constant = 2 #A
        v = 0.4 # Exponent derived empirically Ref. 1 for monomers and oligomers
        Rg = R_constant * num_atoms**v
        
        
    else:
        print("Input error!")
        return 0
    
    
    print("Number of atoms: {} \nRadius of Gyration: {:.2f}".format(num_atoms,Rg))
    if method.lower() == 'rosenthal_henderson':
        d_cutoff = 2*np.pi*Rg
        f_cutoff = 1/d_cutoff
    elif method.lower() == 'singer':
        ko = num_atoms**(-1/12) # Guiner transition non-dimensional
        
        Ro = Rg * np.cbrt(1/num_atoms) # Unit cell dimension around each atom
        
        f_cutoff = ko/Ro
        d_cutoff = 1/f_cutoff
    
    print("Frequency cutoff: {:.2f} (in 1/A) \n".format(f_cutoff))
    print("Frequency cutoff: {:.2f} (in A) \n ".format(d_cutoff))
    
    if return_as_frequency:
        return f_cutoff
    else:
        return d_cutoff



def get_atomic_positions_between_residues(gemmi_structure, chain_name, res_range = None):
    '''
    Extract atom positions between residue range

    Parameters
    ----------
    gemmi_structure : gemmi.Structure()
        input gemmi structure
    chain_name : str
        
    res_range : list
        res_range = [start_res, end_res] (both incl)

    Returns
    -------
    pdb_positions : list
    
    pdb_positions = [[x1, y1, z1], [x2, y2, z3]...] (values in Angstorm)
    '''
    gemmi_model = gemmi_structure[0]

    pdb_positions = []
    for chain in gemmi_model:
        if chain.name == chain_name:
            if res_range is not None:
                for res in chain:
                    if res.seqid.num >= res_range[0] and res.seqid.num <= res_range[1] :
                        for atom in res:
                            pdb_positions.append(atom.pos.tolist())
            else:
                for res in chain:
                    for atom in res:
                        pdb_positions.append(atom.pos.tolist())
                        
    
    return pdb_positions
