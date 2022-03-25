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

def replace_pdb_column_with_arrays(input_pdb, replace_column, replace_array):
    '''
    Replace a column in the PDB (either bfactor or occupancy) with values from an array where the array index position correspods to atom location in the cra generator

    Parameters
    ----------
    input_pdb : TYPE
        DESCRIPTION.
    replace_column : TYPE
        DESCRIPTION.
    replace_array : TYPE
        DESCRIPTION.

    Returns
    -------
    None.

    '''
    import gemmi
    from locscale.include.emmer.pdb.pdb_to_map import detect_pdb_input
    
    st = detect_pdb_input(input_pdb)
    
    replace_array = np.clip(replace_array, 0, 999)  ## PDB array only allows a maximum of three digits
    for i,cra_gen in enumerate(st[0].all()):
        if replace_column=="bfactor":
            cra_gen.atom.b_iso = replace_array[i]
        elif replace_column=="occ":
            cra_gen.atom.occ = replace_array[i]
        else:
            raise UserWarning("Please input either bfactor or occ for the replace_column variable")
        
    
    return st

        
    

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
    
    
    polar_vector: (r, theta, phi) !!
    Theta is the angle in the XY plane

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
            r, theta, phi = vector
            x = r * np.sin(phi) * np.cos(theta)
            y = r * np.sin(phi) * np.sin(theta)
            z = r * np.cos(phi)
            cartesians.append(np.array([x,y,z]))
        return np.array(cartesians)
    else:
        r, theta, phi = polar_vector
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

def get_random_polar_vector(rmsd_magnitude, randomisation="uniform", mean=None):
    if randomisation == "normal":
        if mean is not None:
            r = abs(np.random.normal(loc=mean, scale=rmsd_magnitude))  ## r will be a normally distributed, positive definite variable
            theta = np.random.uniform(0, 2*np.pi)
            phi = np.random.uniform(0, np.pi)
        else:
            r = abs(np.random.normal(loc=0, scale=rmsd_magnitude))  ## r will be a normally distributed, positive definite variable
            theta = np.random.uniform(0, 2*np.pi)
            phi = np.random.uniform(0, np.pi)
                                        
    elif randomisation == "uniform":
        r = np.random.uniform(low=0, high=rmsd_magnitude*2)
        theta = np.random.uniform(0, 2*np.pi)
        phi = np.random.uniform(0, np.pi)
    else:
        raise ValueError("The variable randomisation has only two inputs: normal or uniform")
    
    return np.array([r, theta, phi])

def check_position_inside_mask(position, mask_data):
    value_at_position = mask_data[position[2],position[1], position[0]]
    if value_at_position > 0.9:
        return True
    else:
        return False

def compute_rmsd_two_pdb(input_pdb_1, input_pdb_2, use_gemmi_structure=True, return_array=False):
    from locscale.include.emmer.pdb.pdb_to_map import detect_pdb_input
    from scipy.spatial import distance
    from locscale.include.emmer.pdb.pdb_tools import get_all_atomic_positions
    
    if use_gemmi_structure:
        st_1 = detect_pdb_input(input_pdb_1)
        st_2 = detect_pdb_input(input_pdb_2)
        
        num_atoms_1 = st_1[0].count_atom_sites()
        num_atoms_2 = st_2[0].count_atom_sites()
        
        positions_1 = get_all_atomic_positions(st_1)
        positions_2 = get_all_atomic_positions(st_2)
        
        assert num_atoms_1 == num_atoms_2
    
    else:
        positions_1 = input_pdb_1
        positions_2 = input_pdb_2
        
        assert len(positions_1) == len(positions_2)
        
    
    atomic_distance = []
    for index in range(len(positions_1)):
        dist = distance.euclidean(positions_1[index], positions_2[index])
        atomic_distance.append(dist)
    
    atomic_distance = np.array(atomic_distance)
    
    if return_array:
        return atomic_distance
    else:
        return np.mean(atomic_distance)

def check_mrc_indexing(input_mask, threshold=0.9):
    from locscale.include.emmer.ndimage.map_utils import parse_input
    from locscale.include.emmer.ndimage.map_utils import get_all_voxels_inside_mask

    mask_data = parse_input(input_mask)
    
    binary_mask = (mask_data>=threshold).astype(np.int_)
    
    zero_data = np.zeros(mask_data.shape)
    voxels_inside_mask = set([tuple(x) for x in get_all_voxels_inside_mask(mask_input=mask_data, mask_threshold=threshold)])  ## ZYX format
    
    for voxel in voxels_inside_mask:
        zero_data[voxel[0],voxel[1],voxel[2]] = 1
    
    new_binary_mask = zero_data - binary_mask
    
    assert new_binary_mask.any() == 0 
    
    print("================== OK =========================")
    
    
    
def test_coordinate_functions(rmsd_magnitude=15):
    from locscale.include.emmer.pdb.pdb_utils import convert_polar_to_cartesian, convert_cartesian_to_polar
    from locscale.include.emmer.ndimage.map_utils import convert_mrc_to_pdb_position, convert_pdb_to_mrc_position
    from scipy.spatial import distance
    
    for trial in range(1000):
        pos = np.random.uniform(0,500,size=3)
        r = np.random.uniform(0, rmsd_magnitude*2)
        shake_thetas = np.random.uniform(0, 2*np.pi)
        shake_phis = np.random.uniform(0, np.pi)
        shake_vector = np.column_stack((r, shake_thetas, shake_phis))
        shake_v = shake_vector[0]
        new_pos = pos+convert_polar_to_cartesian(shake_v)
        d = distance.euclidean(new_pos, pos)
        if abs(d-r) > 0.01:
            raise UserWarning("Problem with convert polar to cartesian coordinate")
        
    
    
    for i in range(1000):
        apix = np.random.uniform(0.1,2)
        pos = np.random.uniform(0,500,size=3)
        mrcpos = convert_pdb_to_mrc_position([pos],apix)[0]
        pdbpos = convert_mrc_to_pdb_position([mrcpos],apix)[0]
        d = distance.euclidean(pdbpos, pos)
        if d > apix:
            raise UserWarning("Problem with converting mrc to pdb positions! Getting maximum distance: {} when apix was {}".format(round(d,2), round(apix,2)))
    
    
    print("================== OK =========================")
    
def get_atomic_point_map(mrc_positions, mask_shape):
    zero_map = np.zeros(mask_shape)
    for mrc in mrc_positions:
        zero_map[mrc[0],mrc[1],mrc[2]] = 1
    
    return zero_map

def pick_random_point_within_sphere_of_influence(center_atom, binarised_mask_full, radius_of_influence, apix):
    from locscale.include.emmer.ndimage.map_utils import convert_pdb_to_mrc_position, convert_mrc_to_pdb_position, dilate_mask
    import random
    
    mrc_position = convert_pdb_to_mrc_position([center_atom], apix)[0]
    single_atomic_point_map = get_atomic_point_map([mrc_position], binarised_mask_full.shape)
    
    radius_of_influence_int = int(round(radius_of_influence/apix))
    sphere_of_influence = dilate_mask(binarised_mask_full, radius_of_influence_int, iterations=1)
    
    ## ensure only 0 and 1 exists in the masks
    
    sphere_of_influence = (sphere_of_influence==1).astype(np.int_)
    binarised_mask_full = (binarised_mask_full==1).astype(np.int_)
    
    ## Get intersection space
    
    space_of_influece = sphere_of_influence * binarised_mask_full
    random_mrc_position = random.choice(np.argwhere(space_of_influece==1))
    
    random_pdb_position = convert_mrc_to_pdb_position([random_mrc_position], apix)[0]
    
    return random_pdb_position

def pick_random_point_within_range_kdtree(center_atom, list_of_points_in_mask, range_distance):
    '''
    Pick one point within a distance from a set of points

    Parameters
    ----------
    center_atom : TYPE
        DESCRIPTION.
    list_of_points_in_mask : TYPE
        DESCRIPTION.
    range_distance : TYPE
        DESCRIPTION.

    Returns
    -------
    None.

    '''
    from sklearn.neighbors import KDTree
    tree = KDTree(list_of_points_in_mask)
    
    
    
    
    
def shake_pdb_within_mask(pdb_path, mask_path, rmsd_magnitude, use_pdb_mask=True, masking="strict", threshold=0.5):
    from locscale.include.emmer.pdb.pdb_to_map import detect_pdb_input
    from locscale.include.emmer.ndimage.map_utils import parse_input
    from locscale.include.emmer.pdb.pdb_tools import get_all_atomic_positions, set_all_atomic_positions
    from locscale.include.emmer.ndimage.map_utils import convert_pdb_to_mrc_position, get_all_voxels_inside_mask, convert_mrc_to_pdb_position
    from locscale.include.emmer.ndimage.map_tools import get_atomic_model_mask
    from tqdm import tqdm
    import mrcfile
    import random
    
    ## Inputs
    st = detect_pdb_input(pdb_path)
    
    # If required, use PDB mask with the same shape and apix as a supplied mask_path
    if use_pdb_mask:
        mask = parse_input(get_atomic_model_mask(emmap_path=mask_path, pdb_path=pdb_path, dilation_radius=3))
    else:
        mask = parse_input(mask_path)
    
    outside_mask = np.logical_not(mask>=threshold)
    apix = mrcfile.open(mask_path).voxel_size.tolist()[0]
    
    
    voxel_positions_mask = get_all_voxels_inside_mask(mask, mask_threshold=0.5)
    pdb_positions_mask = convert_mrc_to_pdb_position(voxel_positions_mask, apix)
    
       
    ## Get all atomic positions, in real units (XYZ)
    atomic_positions_values = get_all_atomic_positions(st)  ## XYZ
    # Get shake vector     
    print("Shaking the input structure with input RMSD of {}...".format(rmsd_magnitude))
    shake_radii = np.random.uniform(0, rmsd_magnitude*2, size=len(atomic_positions_values))
    shake_thetas = np.random.uniform(0, 2*np.pi, size=len(atomic_positions_values))
    shake_phis = np.random.uniform(0, np.pi, size=len(atomic_positions_values))
    
    shake_vectors_polar = np.column_stack((shake_radii, shake_thetas, shake_phis))
    print("Shake radii",shake_radii.mean())
    # Get shaken positions in both real units (XYZ) and MRC units (ZYX)
    print("Checking for atoms outside the mask...")
    shaken_atomic_position_native = list(atomic_positions_values + convert_polar_to_cartesian(shake_vectors_polar, multiple=True))
    rmsd_native = compute_rmsd_two_pdb(shaken_atomic_position_native, atomic_positions_values, use_gemmi_structure=False)
    
    shaken_atomic_position= shaken_atomic_position_native.copy()
    
    shaken_mrc_position_list = [tuple(x) for x in convert_pdb_to_mrc_position(shaken_atomic_position, apix)]  ## ZYX
    
    binarised_mask = (mask>=threshold).astype(np.int_)
    mrc_point_map = get_atomic_point_map(shaken_mrc_position_list, binarised_mask.shape)
    available_voxels = get_all_voxels_inside_mask(binarised_mask - mrc_point_map, 1)  
       
    
    random_sample_mrc = random.sample(available_voxels, len(shaken_mrc_position_list))
    random_sample_pdb = convert_mrc_to_pdb_position(random_sample_mrc, apix)
    
    ## Find MRC positions outside the mask
#    mrc_positions_inside_mask = shaken_mrc_position.intersection(voxels_inside_mask)
    np_array_mask_pdb = np.array(pdb_positions_mask)
    
    
    num_atoms_outside = 0
    if masking == "strict":
        from sklearn.neighbors import KDTree
        tree = KDTree(np_array_mask_pdb)
    
    for i,mrc_pos in enumerate(tqdm(shaken_mrc_position_list, "Validating positions")):
        if outside_mask[mrc_pos[0],mrc_pos[1],mrc_pos[2]]:
            if masking == "strict":
                
                #random_point_in_available_space = pick_random_point_within_sphere_of_influence(center_atom=atomic_positions_values[i],
                #                                                                                   binarised_mask_full=binarised_mask,
                #                                                                                   radius_of_influence=rmsd_magnitude*2,
                #                                                                                   apix=apix)
                
                
                neighborhood_indices_list = tree.query_radius(shaken_atomic_position[i:i+1], r=rmsd_magnitude*2)
                
                random_index = random.choice(list(neighborhood_indices_list)[0])
                random_position = np_array_mask_pdb[random_index] + np.random.uniform(0,apix/2,3)
                
                shaken_atomic_position[i] = random_position
            else:
                shaken_atomic_position[i] = random_sample_pdb[i]  
            
            num_atoms_outside += 1
        
    print("{} atoms found outside the mask! Randomly placing them inside the mask:".format(num_atoms_outside))
           
        ## For every mrc position in the like mrc_positions_outside_mask, find the corresponding index and replace it with a value from random_voxels
    ## Convert the list into a dictionary
    
    print("Done... Now converting into PDB")
    shaken_atomic_positions_dictionary = {}
    for i,atomic_position in enumerate(shaken_atomic_position):
        shaken_atomic_positions_dictionary[i] = atomic_position
    
    
    shaken_structure = set_all_atomic_positions(st, shaken_atomic_positions_dictionary)
    
    

    
    rmsd = compute_rmsd_two_pdb(st, shaken_structure)
    print("RMSD between input structure and native shaken structure: {} A".format(round(rmsd_native,2)))
    print("RMSD between the input and output structure is: {} A".format(round(rmsd,2)))
    
    return shaken_structure

def test_pdb_within_mask(input_pdb, mask_path, mask_threshold):
    from locscale.include.emmer.pdb.pdb_to_map import detect_pdb_input
    from locscale.include.emmer.ndimage.map_utils import convert_pdb_to_mrc_position, get_all_voxels_inside_mask, convert_mrc_to_pdb_position
    from locscale.include.emmer.pdb.pdb_tools import get_all_atomic_positions
    import mrcfile

    st = detect_pdb_input(input_pdb)
    mask_data = mrcfile.open(mask_path).data
    apix = mrcfile.open(mask_path).voxel_size.tolist()[0]
        
    binarised_mask = (mask_data>=mask_threshold).astype(np.int_)
    
    pdb_positions = get_all_atomic_positions(st)
    mrc_positions = convert_pdb_to_mrc_position(pdb_positions, apix)
    
    count = 0
    for mrc in mrc_positions:
        count += binarised_mask[mrc[0],mrc[1],mrc[2]]
    
    num_atoms = len(pdb_positions)
    
    percentage_capture = count / num_atoms * 100
    
    print("Percentage capture of PDB within mask at given threshold is {}".format(round(percentage_capture,2)))
    
    return percentage_capture

    
    

    
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
 
                    shake_vector_polar = get_random_polar_vector(rmsd_magnitude=magnitude, randomisation=randomisation)
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


