# -*- coding: utf-8 -*-
"""
Created on Sun Jul 18 14:58:49 2021

"""

# convert_to_mmcif contains a function to convert .pdb file format to .mmcif
# file format using gemmi, as well as functions for the saving of a PDB 
# structure in either format

# global imports
import gemmi

#%% functions
def convert_pdb_to_mmcif(pdb_path, mmcif_path=None):
    """Convert coordinate file from PDB format to mmCIF/PDBx format"""
    # local imports
    import os
    
    structure = gemmi.read_structure(pdb_path)
    if mmcif_path is None:
        mmcif_path = os.path.splitext(pdb_path)[0] + '.cif'
    write_structure_as_mmcif(structure, mmcif_path)
    return mmcif_path


def write_structure_as_mmcif(structure, mmcif_name):
    """Write a Gemmi structure out to an mmCIF file."""
    structure.make_mmcif_document().write_file(mmcif_name)


def write_structure_as_pdb(structure, pdb_name):
    """Write a Gemmi structure out to a PDB file."""
    structure.write_pdb(pdb_name)



