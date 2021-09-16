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


