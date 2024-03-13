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
        gemmi_structure.make_mmcif_document().write_file(pdb_path)
    else:
        gemmi_structure.make_mmcif_document().write_file(new_pdb_path)
        
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
        gemmi_st = gemmi.read_structure(model_path)
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
        gemmi_st = gemmi.read_structure(model_path)
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
        num_residues = num_atoms/8
        Rg = R_constant * num_residues**v
        protein_density = 0.8 ## 0.8 dalton/ang^3 from Henderson, 1995
        molecular_weight = mask_vol_A3 * protein_density
    elif mask is not None and apix is not None and mask_path is None:
        mask_vol_A3, protein_mass, num_atoms, mask_dims,maskshape = measure_mask_parameters(mask=mask, apix=apix, detailed_report=True, verbose=False)
        
        R_constant = 2 #A
        v = 0.4 # Exponent derived empirically Ref. 1 for monomers and oligomers
        num_residues = num_atoms/8
        Rg = R_constant * num_residues**v
        protein_density = 0.8 ## 0.8 dalton/ang^3 from Henderson, 1995
        molecular_weight = mask_vol_A3 * protein_density
    elif num_atoms is not None:
         mol_weight = num_atoms * 13  # daltons 
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


def combine_pdb_structures_into_one(list_of_input_pdb, return_chain_counts=False):
    import gemmi
    from locscale.include.emmer.pdb.pdb_to_map import detect_pdb_input

    def add_chains(combined_model, model_to_add, starting_chain_num):
        import string
        chain_letters_1 = list(string.ascii_uppercase) + list(string.ascii_lowercase)
        # for longer chains, make chain letters as AA, AB, AC, etc.
        chain_letters_2 = [x+y for x in chain_letters_1 for y in chain_letters_1]
        chain_letters = chain_letters_1 + chain_letters_2
        chain_count = starting_chain_num
        for chain in model_to_add:
            combined_model.add_chain(chain_letters[chain_count])
            res_count = 0
            for res in chain:
                atom_count = 0
                residue_name = res.name
                combined_model = add_residue(combined_model, chain_count, res_count, residue_name)
                unique_atoms = get_unique_atoms_in_residue(res)
                for atom in unique_atoms:
                    atom_name = atom.name
                    atom_position = atom.pos
                    atom_element = atom.element
                    combined_model = add_atom(combined_model, chain_count, res_count, atom_count, atom_position, atom_element, atom_name)
                    atom_count += 1

                res_count += 1
            chain_count += 1
        
        final_chain_count = chain_count
        return combined_model, final_chain_count

    def add_residue(gemmi_model, chain_num, res_num, res_name):
        gemmi_model[chain_num].add_residue(gemmi.Residue(),res_num)
        gemmi_model[chain_num][res_num].name = res_name
        gemmi_model[chain_num][res_num].seqid.num = res_num

        return gemmi_model

    def add_atom(gemmi_model, chain_num, res_num, atom_num, gemmi_atom_position,  gemmi_element, atom_name):
        
        atom = gemmi.Atom()
        atom.pos = gemmi_atom_position
        atom.element = gemmi_element
        atom.b_iso = 40
        atom.occ = 1
        atom.name = atom_name

        gemmi_model[chain_num][res_num].add_atom(atom, atom_num)
        
        return gemmi_model
    
    def get_unique_atoms_in_residue(res):
        atoms = []
        atom_names = []
        for atom in res:
            if atom.name not in atom_names:
                atoms.append(atom)
                atom_names.append(atom.name)
        
        return atoms

    
    combined_structure = gemmi.Structure()
    combined_model = gemmi.Model("combined")
    
    chain_count = 0
    final_chain_counts = []
    for input_pdb in list_of_input_pdb:
        chain_count_init = chain_count
        st = detect_pdb_input(input_pdb)
        model_to_add = st[0]
        combined_model, chain_count_final = add_chains(combined_model, model_to_add, chain_count_init)
        chain_count = chain_count_final
        final_chain_counts.append(chain_count_final)

    combined_structure.add_model(combined_model)

    total_num_atoms_in_combined_structure = combined_structure[0].count_atom_sites()

    num_atoms_in_input = 0
    for input_pdb in list_of_input_pdb:
        st = detect_pdb_input(input_pdb)
        model_to_add = st[0]
        num_atoms_in_input += model_to_add.count_atom_sites()
    
    assert total_num_atoms_in_combined_structure == num_atoms_in_input, \
        "Number of atoms in combined structure does not match number of atoms in input"
    
    if return_chain_counts:
        return combined_structure, final_chain_counts
    else:
        return combined_structure

def add_pseudoatoms_to_input_pdb(pdb_path, mask_path, emmap_path, mask_threshold = 0.5, averaging_window=3, pseudomodel_method = "gradient", pseudomodel_iteration=50, bond_length=1.2, fsc_resolution=None, return_chain_counts=False, return_difference_mask=False):
    from locscale.preprocessing.headers import run_pam
    from locscale.include.emmer.ndimage.map_utils import measure_mask_parameters, load_map, save_as_mrc
    from locscale.include.emmer.pdb.pdb_tools import combine_pdb_structures_into_one
    from locscale.include.emmer.ndimage.map_tools import find_unmodelled_mask_region
    import os

    mask, apix = load_map(mask_path)
    # Get the difference mask 
    difference_mask = find_unmodelled_mask_region(fdr_mask_path = mask_path, pdb_path = pdb_path, fdr_threshold = mask_threshold, \
        atomic_mask_threshold = 0.5, averaging_window_size = averaging_window, fsc_resolution = fsc_resolution)
    
    difference_mask_path_filename = mask_path[:-4] + "_difference_mask.mrc"
    difference_mask_path = os.path.join(os.path.dirname(mask_path), difference_mask_path_filename)
    save_as_mrc(difference_mask, difference_mask_path, apix)

    num_atoms, _ = measure_mask_parameters(mask_path = difference_mask_path, edge_threshold=mask_threshold, verbose=False)
    if num_atoms == 0:
        raise ValueError("No unmodelled atoms found in difference mask. Run Model Based LocScale without passing the '--complete_model' flag.")
    
    pdb_filename = os.path.basename(pdb_path)
    print("Adding {} pseudoatoms to {}".format(num_atoms, pdb_filename))

    partial_pseudo_model_path = run_pam(emmap_path = emmap_path, mask_path = difference_mask_path, threshold = mask_threshold, \
       num_atoms = num_atoms, method=pseudomodel_method, bl=bond_length, total_iterations = pseudomodel_iteration, verbose=False)

    combined_structure, final_chain_counts = combine_pdb_structures_into_one([pdb_path, partial_pseudo_model_path], return_chain_counts=return_chain_counts)

    return_elements = [combined_structure]
    if return_chain_counts:
        return_elements.append(final_chain_counts)
    if return_difference_mask:
        return_elements.append(difference_mask_path)
    
    if not return_chain_counts and not return_difference_mask:
        return combined_structure
    else:
        return tuple(return_elements)
def get_coordinate_bfactors(st):
    atom_info = []
    for cra in st[0].all():
        atom = cra.atom
        atom_info.append((atom.pos, atom.b_iso))
    
    return atom_info

def check_if_atomic_position_in_mask(atomic_pos, mask, apix):
    x, y, z = atomic_pos
    int_x, int_y, int_z = int(round(x/apix)), int(round(y/apix)), int(round(z/apix))
    return mask[int_z, int_y, int_x]


def neighborhood_bfactor_correlation_sample(input_pdb, min_radius=1, max_radius=10, num_steps=10, sample_size=1000, mask_path=None, invert=False):
    import gemmi
    from locscale.include.emmer.pdb.pdb_to_map import detect_pdb_input
    from locscale.include.emmer.pdb.pdb_utils import get_coordinates, get_bfactors
    from locscale.include.emmer.ndimage.map_utils import convert_pdb_to_mrc_position, load_map
    from tqdm import tqdm
    import random
    from scipy.stats import pearsonr
    st = detect_pdb_input(input_pdb)
    atom_info = get_coordinate_bfactors(st)
    if mask_path is not None:
        mask, apix = load_map(mask_path)
        mask = (mask > 0.5).astype(bool)   
        if invert:
            mask = np.logical_not(mask) 
        atom_info_filtered = [atom for atom in atom_info if check_if_atomic_position_in_mask(atom[0], mask, apix)] 
    else:
        atom_info_filtered = atom_info

    sample_atom_info = random.sample(atom_info_filtered, sample_size)

    bfactor_correlation_with_distance = {}
    for r in tqdm(np.linspace(min_radius,max_radius,num_steps), total=num_steps, desc="Finding neighborhood bfactor correlation"):
        sample_bfactors = []
        sample_coordinates = []
        ns_pseudo = gemmi.NeighborSearch(st[0], st.cell, r).populate()
        neighborhood_bfactor_list = []
        for sample_info in sample_atom_info:
            sample_bfactors.append(sample_info[1])
            sample_coordinates.append(sample_info[0].tolist())
            gemmi_coord = sample_info[0]
            neighbors = ns_pseudo.find_atoms(gemmi_coord, '\0', radius=r)
            if len(neighbors) > 0:
                neighborhood_bfactor_list.append(np.mean([n.to_cra(st[0]).atom.b_iso for n in neighbors]))
                
        bfactor_correlation_at_r = pearsonr(sample_bfactors, neighborhood_bfactor_list)
        bfactor_correlation_with_distance[r] = [sample_bfactors, neighborhood_bfactor_list, bfactor_correlation_at_r, sample_coordinates]

    return bfactor_correlation_with_distance

def neighborhood_bfactor_correlation(input_pdb, min_radius=1, max_radius=10, num_steps=10, mask_path=None, invert=False):
    import gemmi
    from locscale.include.emmer.pdb.pdb_to_map import detect_pdb_input
    from locscale.include.emmer.ndimage.map_utils import convert_pdb_to_mrc_position, load_map
    from tqdm import tqdm
    from scipy.stats import pearsonr
    st = detect_pdb_input(input_pdb)
    model = st[0]
    if mask_path is not None:
        mask, apix = load_map(mask_path)
        mask = (mask > 0.5).astype(bool)
        if invert:
            mask = np.logical_not(mask)


    bfactor_correlation_with_distance = {}
    for r in tqdm(np.linspace(min_radius,max_radius,num_steps), total=num_steps, desc="Finding neighborhood bfactor correlation"):
        ns_pseudo = gemmi.NeighborSearch(st[0], st.cell, r).populate()
        individual_bfactor_list = []
        neighborhood_bfactor_list = []
        atomic_position_list = []
        for cra in model.all():
            atom = cra.atom
            if mask_path is not None:
                if not check_if_atomic_position_in_mask(atom.pos, mask, apix):
                    continue
            atomic_biso = atom.b_iso
            neighbors = ns_pseudo.find_atoms(atom.pos, '\0', radius=r)
            neigbor_atoms = [x.to_cra(model).atom for x in neighbors]
            atomic_bfactor_list = np.array([x.b_iso for x in neigbor_atoms])
            average_bfactor_neighbors = atomic_bfactor_list.mean()
            individual_bfactor_list.append(atomic_biso)
            neighborhood_bfactor_list.append(average_bfactor_neighbors)
            atomic_position_list.append(atom.pos.tolist())
            
        
        pearson_correlation = pearsonr(individual_bfactor_list, neighborhood_bfactor_list)
        bfactor_correlation_with_distance[r] = [individual_bfactor_list,neighborhood_bfactor_list, pearson_correlation, atomic_position_list]
    
    return bfactor_correlation_with_distance
