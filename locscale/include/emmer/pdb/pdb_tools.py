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
#from emmer.ndimage.filter import *
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
    from locscale.include.emmer.pdb.pdb_utils import shift_coordinates
    
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

def find_wilson_cutoff(model_path=None, input_gemmi_st=None, mask_path=None, mask = None, apix=None, num_atoms=None, method='Singer', return_as_frequency=False, verbose=True):
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
    
    from locscale.include.emmer.ndimage.map_utils import measure_mask_parameters
    if model_path is not None:
        gemmi_st = gemmi.read_pdb(model_path)
        num_atoms = gemmi_st[0].count_atom_sites()
        Rg = find_radius_of_gyration(input_gemmi_st=gemmi_st)
        molecular_weight = gemmi_st[0].calculate_mass()
    elif input_gemmi_st is not None:
        gemmi_st = input_gemmi_st.clone()
        num_atoms = gemmi_st[0].count_atom_sites()
        Rg = find_radius_of_gyration(input_gemmi_st=gemmi_st)
        molecular_weight = gemmi_st[0].calculate_mass()
    elif mask_path is not None:
        mask_vol_A3, protein_mass, num_atoms, mask_dims,maskshape = measure_mask_parameters(mask_path=mask_path, detailed_report=True, verbose=False)
        R_constant = 2 #A
        v = 0.4 # Exponent derived empirically Ref. 1 for monomers and oligomers
        Rg = R_constant * num_atoms**v
        protein_density = 0.8 ## 0.8 dalton/ang^3 from Henderson, 1995
        molecular_weight = mask_vol_A3 * protein_density
    elif mask is not None and apix is not None and mask_path is None:
        mask_vol_A3, protein_mass, num_atoms, mask_dims,maskshape = measure_mask_parameters(mask=mask, apix=apix, detailed_report=True, verbose=False)
        
        R_constant = 2 #A
        v = 0.4 # Exponent derived empirically Ref. 1 for monomers and oligomers
        Rg = R_constant * num_atoms**v
        protein_density = 0.8 ## 0.8 dalton/ang^3 from Henderson, 1995
        molecular_weight = mask_vol_A3 * protein_density
    elif num_atoms is not None:
         mol_weight = num_atoms * 16  # daltons 
         wilson_cutoff_local = 1/(0.309 * np.power(mol_weight, -1/12))   ## From Amit Singer
         return wilson_cutoff_local
    
    else:
        print("Input error!")
        return 0
    
    if verbose:
        print("Number of atoms: {} \nRadius of Gyration: {:.2f}\n".format(num_atoms,Rg))
        print("Molecular weight estimated to be {} kDa\n".format(round(molecular_weight/1000,1)))
        
    if method.lower() == 'rosenthal_henderson':
        d_cutoff = 2*np.pi*Rg
        f_cutoff = 1/d_cutoff
    elif method.lower() == 'singer':
        
        f_cutoff = 0.309 * np.power(molecular_weight, -1/12)  ## From Singer, 2021
        d_cutoff = 1/f_cutoff
    
    if verbose:
        print("Frequency cutoff: {:.2f}  (= {:.2f} A resolution)\n".format(f_cutoff, d_cutoff))
    #print("Frequency cutoff: {:.2f} (in A) \n ".format(d_cutoff))
    
    if return_as_frequency:
        return f_cutoff
    else:
        return d_cutoff



def get_all_atomic_positions(gemmi_structure, as_dictionary=False):
    '''
    Extract atom positions

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
    import gemmi
    
    st = gemmi_structure.clone()
    
    if as_dictionary:
        pdb_positions = {}
        for i,cra_obj in enumerate(st[0].all()):
            pdb_positions[i] = np.array(cra_obj.atom.pos.tolist())
        
        return pdb_positions
                        
    
    else:
        pdb_positions = []
        for chain in st[0]:
            for res in chain:
                for atom in res:
                    pdb_positions.append(atom.pos.tolist())
                            
        
        return np.array(pdb_positions)

def set_all_atomic_positions(gemmi_structure, positions_dictionary):
    '''
    Input a dictionary where keys are atomic "access IDs " generated by the function get_all_atomic_positions

    Parameters
    ----------
    gemmi_structure : TYPE
        DESCRIPTION.
    positions_dictionary : TYPE
        DESCRIPTION.

    Returns
    -------
    None.

    '''
    import gemmi
    
    st = gemmi_structure.clone()
    for i, cra_obj in enumerate(st[0].all()):
        new_position = gemmi.Position(positions_dictionary[i][0],positions_dictionary[i][1],positions_dictionary[i][2])
        cra_obj.atom.pos = new_position
    
    return st
    


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

def find_number_of_neighbors(input_pdb, atomic_position, window_size_A):
    from locscale.include.emmer.pdb.pdb_to_map import detect_pdb_input
    import gemmi
    
    
    gemmi_structure = detect_pdb_input(input_pdb)
    
    # Neighbor Search initialize
    
    ns = gemmi.NeighborSearch(gemmi_structure[0], gemmi_structure.cell, window_size_A).populate()
    
    gemmi_position = gemmi.Position(atomic_position[0], atomic_position[1], atomic_position[2])
    
    neighbors = ns.find_atoms(gemmi_position, '\0', radius=window_size_A)
    atoms = [gemmi_structure[0][x.chain_idx][x.residue_idx][x.atom_idx] for x in neighbors]
    number_of_neighbors = len(atoms)
    
    
    return number_of_neighbors

def get_atomic_bfactor_window(input_pdb, atomic_position, window_size_A, min_dist=0.1):
    from locscale.include.emmer.pdb.pdb_to_map import detect_pdb_input
    import gemmi
    
    
    gemmi_structure = detect_pdb_input(input_pdb)
    
    # Neighbor Search initialize
    
    ns = gemmi.NeighborSearch(gemmi_structure[0], gemmi_structure.cell, window_size_A).populate()
    
    gemmi_position = gemmi.Position(atomic_position[0], atomic_position[1], atomic_position[2])
    
    neighbors = ns.find_atoms(gemmi_position, '\0', radius=window_size_A)
    atoms = [gemmi_structure[0][x.chain_idx][x.residue_idx][x.atom_idx] for x in neighbors]
    atomic_bfactor_list = np.array([x.b_iso for x in atoms])
    
    average_atomic_bfactor = atomic_bfactor_list.mean()
    
    return average_atomic_bfactor


    
    