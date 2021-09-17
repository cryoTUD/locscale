# -*- coding: utf-8 -*-
"""
Created on Sun Jul 18 14:52:11 2021
"""

# get_dssp contains several functions to analyse secondary structure
# information in PDB structures using Biopython

# global imports
from Bio.PDB import PDBParser
from Bio.PDB.DSSP import DSSP
from Bio.PDB import PDBList

#%% functions

def get_dssp(pdbid,pdb_path):     
    dsspfile_location = '/home/alok/soft/ccp4-7.1/ccp4-7.1/bin/mkdssp'
    parser = PDBParser()
    structure = parser.get_structure(pdbid,pdb_path)
    model = structure[0]
    
    try:
         dssp = DSSP(model,pdb_path,dssp=dsspfile_location)
         print("Successfully performed DSSP! \n")
         
         
    except:
         print("Problem with DSSP")
         return 0
    
    return dssp
 
def get_secondary_structure_residues_content(pdbid,pdb_path=None):
     '''
     Input: pdbid: string '3j5p' 
     pdb_path: string - gives the full path of the PDB file
     '''
     # local imports
     import os

     deletefile = False
     if pdb_path is None:
         deletefile = True
         pdbl = PDBList()
         localfile = pdbl.retrieve_pdb_file(pdbid,file_format='pdb')
         print("Successfully downloaded PDB! \n")
         pdb_path = localfile
     
     dssp = get_dssp(pdbid,pdb_path)    
     secondary_structures = [dssp[key][2] for key in dssp.keys()]
     helix_residues = secondary_structures.count('H') + secondary_structures.count('G') + secondary_structures.count('I')
     sheet_residues = secondary_structures.count('B') + secondary_structures.count('E')
     loop_residues = secondary_structures.count('T') + secondary_structures.count('S')
     total_residues = len(secondary_structures)
     
     secondary_structure_distribution = {}
     secondary_structure_distribution['helix'] = helix_residues/total_residues
     secondary_structure_distribution['sheet'] = sheet_residues/total_residues
     secondary_structure_distribution['loop'] = loop_residues/total_residues
     secondary_structure_distribution['total'] = total_residues
     
     if deletefile:
         os.remove(localfile)
         directory = pdbid[1:3]
         os.rmdir(directory+'/')
     
     return secondary_structures,secondary_structure_distribution
    
def link_gemmi_model_with_dssp(gemmi_model,dssp):
     '''
     gemmi_model: gemmi.Model() 
     dssp: dictionary which is the output of Bio.PDB.DSSP.DSSP() function
     
     return:
     linked_dictionary[(chain_name,res_id)] = relevant_dssp_key
     '''
     linked_dictionary = {}
     for key in dssp.keys():
          chain_name = key[0]
          res_id = key[1][1]
          linked_dictionary[(chain_name,res_id)] = key
          
     return linked_dictionary

def split_gemmi_model_based_on_dssp(pdbid,pdb_path):
     
     '''
     Returns two gemmi.Model() based on secondary structure
     '''
     #local imports
     import gemmi
     
     dssp = get_dssp(pdbid,pdb_path)
     secondary_structure_res_list,_ = get_secondary_structure_residues_content(pdbid,pdb_path)
     
     gemmi_model = gemmi.read_pdb(pdb_path)[0]
     helix_model = gemmi.Model('H')
     sheet_model = gemmi.Model('S')
     skipped_residues = 0
     
     linked_dictionary = link_gemmi_model_with_dssp(gemmi_model,dssp)
     print("Number of dssp keys >>> " + str(len(dssp.keys())))
     
     num_residue = 0
     for chain in gemmi_model:
          helix_model.add_chain(chain.name)
          sheet_model.add_chain(chain.name)

          for res in chain:
               #dssp_key = tuple([chain.name,tuple([' ',res.seqid.num,' '])])
               num_residue += 1
               try:
                    dssp_key = linked_dictionary[(chain.name,res.seqid.num)]
                    if dssp[dssp_key][2] in ['H','G','I']:
                         helix_model[chain.name].add_residue(res)
                    elif dssp[dssp_key][2] in ['B','E']:
                         sheet_model[chain.name].add_residue(res)
               except KeyError:
                    skipped_residues += 1
                    
     print("Number of residues in Gemmi Model  >>> " + str(num_residue))

     return helix_model, sheet_model,tuple([num_residue,skipped_residues])

