# -*- coding: utf-8 -*-
"""
Created on Sun Jul 18 14:08:55 2021
"""

# pdb_util contains helper functions used in other applications in the emmer
# toolbox. pdb_util functions are classified as such if there is little need
# to manually call these functions and/or when their output is used in
# many different applications.

# global imports
import gemmi
import numpy as np

#%% functions
     
def check_if_gemmi_st_is_downloadable(pdb_id):
    """[summary]

    Args:
        pdb_id (str): PDB ID, like: "3j5p"

    Returns:
        Bool: indicating whether the gemmi structure is downloadable (True) or not (False)
    """    
    import pypdb
    
    try: 
        pdb_file = pypdb.get_pdb_file(pdb_id, filetype='pdb',compression=False)
        return True
    except AttributeError:
        try:
            cif_file = pypdb.get_pdb_file(pdb_id, filetype='cif', compression=False)
            return True
        except Exception as e:
            print("- Exception: {}".format(e))
            return False

def get_gemmi_st_from_id(pdb_id):
    '''
    Returns a gemmi.Structure() containing the PDB coordinate information and headers for a given PDB-ID. 

    Parameters
    ----------
    pdb_id : string
        PDB ID: "3j5p"

    Returns
    -------
    gemmi_st : gemmi.Structure()

    '''          
    import pypdb
    import gemmi
    import os
    from gemmi import cif
    
    try: 
        pdb_file = pypdb.get_pdb_file(pdb_id,filetype='pdb',compression=False)
        pdb_struc = gemmi.read_pdb_string(pdb_file)
        print("- Successfully downloaded PDBgemmi {} from database".format(pdb_id))
        return pdb_struc
    except AttributeError:
        cif_file = pypdb.get_pdb_file(pdb_id, filetype='cif',\
                                          compression=False)
        doc_file = cif.read_string(cif_file)
        doc_file.write_file('cif_file.cif')
        cif_struc = gemmi.read_structure("cif_file.cif")
        os.remove("cif_file.cif")
        return cif_struc

def shift_coordinates(in_model_path=None, out_model_path=None,\
                      trans_matrix=[0,0,0], remove_charges=False,\
                      input_structure=None):
    """
    Shift atomic coordinates based on a translation matrix
    """
    if input_structure is None:
        structure = gemmi.read_structure(in_model_path)
    else: 
        structure = input_structure.clone()
        
    for model in structure:
        for chain in model:
            for residue in chain:
                for atom in residue:
                    atom.pos = gemmi.Position(atom.pos.x+trans_matrix[0],
                                              atom.pos.y+trans_matrix[1],
                                              atom.pos.z+trans_matrix[2])
                    if remove_charges: #remove charges
                        if atom.charge != 0:
                            atom.charge = 0
    if input_structure is None:
        structure.write_pdb(out_model_path)
    else:
        return structure
  
def split_model_based_on_nucleotides(gemmi_st):
    dna_model = gemmi.Model('D')
    rna_model = gemmi.Model('R')
    
    for model in gemmi_st:
        for chain in model:
            dna_model.add_chain(chain.name)
            rna_model.add_chain(chain.name)
            for res in chain:
                if res.name in ['C','G','A','U','I']:
                    rna_model[chain.name].add_residue(res)
                elif res.name in ['DC','DG','DA','DU','DI','DT']:
                    dna_model[chain.name].add_residue(res)
    
    dna_st = gemmi.Structure()
    dna_st.add_model(dna_model)
    
    rna_st = gemmi.Structure()
    rna_st.add_model(rna_model)

    return dna_st, rna_st    

def convert_polar_to_cartesian(polar_vector, multiple=False):
    '''
    Convert polar to cartesian.. Blindly following the formula mentioned here: 
        https://math.libretexts.org/Bookshelves/Calculus/Book%3A_Calculus_(OpenStax)/12%3A_Vectors_in_Space/12.7%3A_Cylindrical_and_Spherical_Coordinates#:~:text=To%20convert%20a%20point%20from,and%20z%3D%CF%81cos%CF%86.
        (accessed: 23-2-2022) 
    
    

    Parameters
    ----------
    r : float
        
    phi : float 
        first angle in radians
    theta : float
        second angle in radians

    Returns
    -------
    cartesian : numpy.ndarray [1x3]
        (x,y,z)

    '''
    if multiple:
        cartesians = []
        for vector in polar_vector:
            r, phi, theta = vector
            x = r * np.sin(phi) * np.cos(theta)
            y = r * np.sin(phi) * np.sin(theta)
            z = r * np.cos(phi)
            cartesians.append(np.array([x,y,z]))
        return np.array(cartesians)
    else:
        r, phi, theta = polar_vector
        x = r * np.sin(phi) * np.cos(theta)
        y = r * np.sin(phi) * np.sin(theta)
        z = r * np.cos(phi)
    
        cartesian = np.array([x,y,z])
    
        return cartesian

def convert_cartesian_to_polar(cartesian):
    '''
    Same as above

    Parameters
    ----------
    x : TYPE
        DESCRIPTION.
    y : TYPE
        DESCRIPTION.
    z : TYPE
        DESCRIPTION.

    Returns
    -------
    Polar : numpy.ndarray

    '''
    x, y, z = cartesian
    r = np.sqrt(np.power(x,2)+np.power(y,2)+np.power(z,2))
    theta = np.arctan(y/x)
    phi = np.arccos(z / r)
    
    polar = np.array([r, theta, phi])
    
    return polar

def get_random_polar_vector(magnitude, randomisation="uniform", mean=None):
    if randomisation == "normal":
        if mean is not None:
            r = abs(np.random.normal(loc=mean, scale=magnitude))  ## r will be a normally distributed, positive definite variable
            theta = np.random.uniform(0, np.pi*2)
            phi = np.random.uniform(0, np.pi*2)
        else:
            r = abs(np.random.normal(loc=0, scale=magnitude))  ## r will be a normally distributed, positive definite variable
            theta = np.random.uniform(0, np.pi*2)
            phi = np.random.uniform(0, np.pi*2)
                                        
    elif randomisation == "uniform":
        r = np.random.uniform(low=0, high=magnitude)
        theta = np.random.uniform(0, np.pi*2)
        phi = np.random.uniform(0, np.pi*2)
    else:
        raise ValueError("The variable randomisation has only two inputs: normal or uniform")
    
    return np.array([r, phi, theta])

def check_position_inside_mask(position, mask_data):
    value_at_position = mask_data[position[2],position[1], position[0]]
    if value_at_position > 0.9:
        return True
    else:
        return False
    
def shake_pdb_within_mask(input_pdb, input_mask, magnitude, apix=None):
    from locscale.include.emmer.pdb.pdb_to_map import detect_pdb_input
    from locscale.include.emmer.ndimage.map_utils import parse_input
    from locscale.include.emmer.pdb.pdb_tools import get_all_atomic_positions
    from locscale.include.emmer.ndimage.map_utils import convert_pdb_to_mrc_position, get_all_voxels_inside_mask
    import os
    import mrcfile
    import random
    
    st = detect_pdb_input(input_pdb)
    mask = parse_input(input_mask)
    voxels_inside_mask = set(tuple(get_all_voxels_inside_mask(mask_input=mask, mask_threshold=0.99)))
    if os.path.exists(str(input_mask)):
        apix = mrcfile.open(input_mask).voxel_size.tolist()[0]
    else:
        if apix is None:
            raise UserWarning("Please provide apix")
        else:
            apix = float(apix)
            
    
    atomic_positions = get_all_atomic_positions(st)
    
    
    shake_radii = np.random.uniform(0, magnitude, size=len(atomic_positions))
    shake_phis = np.random.uniform(-np.pi, np.pi, size=len(atomic_positions))
    shake_thetas = np.random.uniform(-np.pi, np.pi, size=len(atomic_positions))
    
    shake_vectors_polar = np.column_stack((shake_radii, shake_phis, shake_thetas))
    
    shaken_atomic_position = atomic_positions + convert_polar_to_cartesian(shake_vectors_polar, multiple=True)
    shaken_mrc_position = set(tuple(convert_pdb_to_mrc_position(shaken_atomic_position, apix)))
    
    atomic_positions_inside_mask = shaken_mrc_position.intersection(voxels_inside_mask)
    
    atomic_positions_outside_mask = list(shaken_mrc_position - voxels_inside_mask)
    new_atomic_positions_to_add = set(tuple(random.sample(voxels_inside_mask, len(atomic_positions_outside_mask))))
    
    final_atomic_positions = list(atomic_positions_inside_mask.union(new_atomic_positions_to_add))
    
    
    
    
    

def shake_pdb(input_pdb, magnitude, randomisation="uniform", mean=None):
    '''
    Function to generate a new pdb by shaking an old PDB

    Parameters
    ----------
    input_pdb : path to pdb or gemmi.Structure()
        DESCRIPTION.
    rmsd : float
        The default is 0.5.

    Returns
    -------
    shaked_st : gemmi.Structure

    '''
    from locscale.include.emmer.pdb.pdb_to_map import detect_pdb_input
    
    input_gemmi_st = detect_pdb_input(input_pdb)
    
    
    assert magnitude > 0
    
    for model in input_gemmi_st:
        for chain in model:
            for res in chain:
                for atom in res:
                    current_pos = np.array(atom.pos.tolist())
 
                    shake_vector_polar = get_random_polar_vector(magnitude=magnitude, randomisation=randomisation)
                    shake_vector_cartesian = convert_polar_to_cartesian(shake_vector_polar)
                        
                    new_pos = current_pos + shake_vector_cartesian
                        
                    atom.pos = gemmi.Position(new_pos[0], new_pos[1], new_pos[2])
    
    return input_gemmi_st

def get_bfactors(in_model_path, return_as_list=True):
    """
    Get B-factors of atoms
    """
    dict_chain = {}
    structure = gemmi.read_structure(in_model_path)
    list_bfact = []
    for model in structure:
        for chain in model:
            polymer = chain.get_polymer()
            #skip non polymers
            #if not polymer: continue
            if not chain.name in dict_chain:
                dict_chain[chain.name] = {}
            for residue in chain:
                
                residue_id = str(residue.seqid.num)+'_'+residue.name
                for atom in residue:
                    list_bfact.append(atom.b_iso)
                avg_bfact = sum(list_bfact)/float(len(list_bfact))
                dict_chain[chain.name][residue_id] = round(avg_bfact,3)
        break # ignore other models
    
    if return_as_list:
        return list_bfact
    else:
        return dict_chain

def add_atomic_bfactors(in_model_path=None, input_gemmi_st=None,
                        additional_biso=None, out_file_path=None):
    '''
    Function to modify atomic bfactors uniformly by adding or subtracting b factors to each atom present in the PDB.
    

    Parameters
    ----------
    in_model_path : str, optional
        Path to a PDB file. 
    gemmi_st : gemmi.Structure()
        Pass a gemmi.Structure() instead of a path to perform computation online
    additional_biso : float, 
        Parameter to specify how the bfactors of the atomic model be modified

    Returns
    -------
    If in_model_path is passed, returns the output model path.. Else returns the gemmi.Structure()

    '''
    
    if in_model_path is not None:
        gemmi_st = gemmi.read_pdb(in_model_path)
    elif input_gemmi_st is not None:
        gemmi_st = input_gemmi_st.clone()
    else:
        print("Pass either the PDB path or the gemmi structure! \n")
        return 0
    
    if additional_biso is not None:
        add_b_iso = additional_biso
    else:
        print("Enter the b factor to add")
        return 0
    
    for model in gemmi_st:
        for chain in model:
            for res in chain:
                for atom in res:
                    atom.b_iso += add_b_iso
    
    if in_model_path is not None:
        if out_file_path is not None:
            output_filepath = out_file_path
        else:
            output_filepath = in_model_path[:-4]+'_modified_bfactor.pdb'
    
        gemmi_st.write_pdb(output_filepath)
    
    else:
        return gemmi_st
    
def set_atomic_bfactors(in_model_path=None, input_gemmi_st=None,
                        b_iso=None, out_file_path=None):
    '''
    Function to modify atomic bfactors uniformly by adding or subtracting b factors to each atom present in the PDB.
    

    Parameters
    ----------
    in_model_path : str, optional
        Path to a PDB file. 
    gemmi_st : gemmi.Structure()
        Pass a gemmi.Structure() instead of a path to perform computation online
    b_iso : float, 
        Parameter to specify the bfactors of the atomic model

    Returns
    -------
    If in_model_path is passed, returns the output model path.. Else returns the gemmi.Structure()

    '''
    
    if in_model_path is not None:
        gemmi_st = gemmi.read_pdb(in_model_path)
    elif input_gemmi_st is not None:
        gemmi_st = input_gemmi_st.clone()
    else:
        print("Pass either the PDB path or the gemmi structure! \n")
        return 0
    
    if b_iso is not None:
        b_iso = b_iso
    else:
        print("Enter the b factor to add")
        return 0
    
    for model in gemmi_st:
        for chain in model:
            for res in chain:
                for atom in res:
                    atom.b_iso = b_iso
    
    if in_model_path is not None:
        if out_file_path is not None:
            output_filepath = out_file_path
        else:
            output_filepath = in_model_path[:-4]+'_modified_bfactor.pdb'
    
        gemmi_st.write_pdb(output_filepath)
    
    else:
        return gemmi_st
    
def calc_bfact_deviation(in_model_path):
    structure = gemmi.read_structure(in_model_path)
    dict_deviation = {}
    for dist in [3.0,5.0,10.0]:
        subcells = gemmi.SubCells(structure[0], structure.cell, dist)
        subcells.populate()
        dict_deviation[dist] = {}
        for model in structure:
            for chain in model:
                polymer = chain.get_polymer()
                #skip non polymers
                if not polymer: continue
                if not chain.name in dict_deviation[dist]:
                    dict_deviation[dist][chain.name] = {}
                for residue in chain.get_polymer():
                    list_bfact = []
                    residue_id = str(residue.seqid.num)+'_'+residue.name
                    for atom in residue:
                        if atom.name == 'CA':
                            ca_bfact = atom.b_iso
                            list_neigh_bfact = []
                            marks = subcells.find_neighbors(atom, min_dist=0.1, max_dist=dist)
                            for mark in marks:
                                cra = mark.to_cra(model)
                                neigh_atom = cra.atom
                                if neigh_atom.name == 'CA':
                                    list_neigh_bfact.append(neigh_atom.b_iso)
                            try: avg_neigh = sum(list_neigh_bfact)/len(list_neigh_bfact)
                            except ZeroDivisionError: pass
                            break
                    if len(list_neigh_bfact) > 0:
                        dict_deviation[dist][chain.name][residue_id] = abs(ca_bfact - avg_neigh)
            break # ignore other models
        
    return dict_deviation

def get_residue_ca_coordinates(in_model_path):
    dict_coord = {}
    structure = gemmi.read_structure(in_model_path)
    for model in structure:
        if not model.name in dict_coord: dict_coord[model.name] = {}
        for chain in model:
            polymer = chain.get_polymer()
            #skip non polymers
            #if not polymer: continue
            if not chain.name in dict_coord[model.name]: 
                dict_coord[model.name][chain.name] = {}
            for residue in chain:
                residue_id = str(residue.seqid.num)+'_'+residue.name
                residue_centre = ()
                if residue.name in ['A','T','C','G','U']:#nuc acid
                    for atom in residue:
                        if atom.name in ["P","C3'","C1'"]:
                            residue_centre = (atom.pos.x,atom.pos.y,atom.pos.z)
                else:
                    for atom in residue:
                        if atom.name == 'CA':#prot
                            residue_centre = (atom.pos.x,atom.pos.y,atom.pos.z)
                if len(residue_centre) == 0:#non nuc acid / prot
                    try: 
                        center_index = len(residue)/2
                        atom = residue[center_index]
                        residue_centre = (atom.pos.x,atom.pos.y,atom.pos.z)
                    except: 
                        for atom in residue:
                            residue_centre = (atom.pos.x,atom.pos.y,atom.pos.z)
                            break #first atom
                if len(residue_centre) > 0:
                    dict_coord[model.name][str(chain.name)][str(residue.seqid.num)] = \
                                            [residue_centre, residue.name]

    return dict_coord

def get_coordinates(in_model_path):
    list_coord = []
    structure = gemmi.read_structure(in_model_path)
    for model in structure:
        for chain in model:
            polymer = chain.get_polymer()
            #skip non polymers
            if not polymer: continue
            for residue in chain:
                residue_id = str(residue.seqid.num)+'_'+residue.name
                for atom in residue:
                    coord = atom.pos #gemmi Position
                    list_coord.append([coord.x,coord.y,coord.z])
    return list_coord

def remove_atomic_charges(in_model_path,out_model_path):
    structure = gemmi.read_structure(in_model_path)
    for model in structure:
        for chain in model:
            polymer = chain.get_polymer()
            #skip non polymers
            #if not polymer: continue
            for residue in chain:
                residue_id = str(residue.seqid.num)+'_'+residue.name
                for atom in residue:
                    if atom.charge != 0:
                        atom.charge = 0
    structure.write_pdb(out_model_path)
    return 1


