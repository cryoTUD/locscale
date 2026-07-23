# -*- coding: utf-8 -*-
"""
Created on Sun Jul 18 15:02:15 2021
"""

# get_associated_pdbs_for_all_emddb contains functions for retrieving all pdb
# id's corresponding to a emdb id. It also allows downloading an density map
# based on a given emdb id.

# global imports
import os
import ftplib
from datetime import datetime

#%% functions

def get_all_pdb_id():
    '''
    Function which returns all PDB ID's from the FTP server
    '''    
    ftp_host = "ftp.wwpdb.org"
    ftp = ftplib.FTP(ftp_host,'*','*')
    directory = "pub/pdb/data/structures/all/pdb/"
    
    ftp.cwd(directory)
    try:
         files = ftp.nlst()
    except:
         print("Oops, error in retreiving file information. Please try again!")
         return 0
    
    pdb_list = [x[3:-7] for x in files]
    ftp.close()
    
    return pdb_list

def get_all_emdb_id():    
    ftp_host = "ftp.ebi.ac.uk"
    ftp = ftplib.FTP(ftp_host,'anonymous')
    directory = "pub/databases/emdb/structures/"
    
    ftp.cwd(directory)
    try:
         files = ftp.nlst()
    except:
         print("Oops, error in retreiving file information. Please try again!")
         return 0
    
    emdb_list = [x[4:] for x in files]
    ftp.close()
    
    return emdb_list

def get_associated_pdbs_for_all_emdb():
    ''' 
    Finds all EMDB filenames. Then for each EMDB file, gets the associated PDB file name.
    Returns a dictionary; associated_pdb['5778'] = ['3j5p','3j5d']
    '''
    # local imports
    from xml.dom import minidom
    
    emdb_list = get_all_emdb_id()
    print("Got all EMDB ID")
    associated_pdbs = {}
    ftp_host = "ftp.ebi.ac.uk"
    ftp = ftplib.FTP(ftp_host,'anonymous')
    directory = "pub/databases/emdb/structures/"
    ftp.cwd(directory)
    for emdb_id in emdb_list:
         print("Now finding PDBs for EMD-"+emdb_id)
    
         header_directory = 'EMD-'+emdb_id+'/header'
         xml_file = 'emd-'+emdb_id+'-v19.xml'
         temp_xml_file = 'temporary.xml'
         ftp.cwd(header_directory)
         ftp.retrbinary("RETR "+xml_file,open(temp_xml_file,'wb').write)
         xmldoc = minidom.parse(temp_xml_file)
         
         pdb_entry_list_xml = xmldoc.getElementsByTagName('fittedPDBEntryIdList')
         pdb_ids = []
         if pdb_entry_list_xml:
              for child in pdb_entry_list_xml[0].childNodes:
                   if child.hasChildNodes():
                        pdb_ids.append(child.childNodes[0].nodeValue)
              associated_pdbs[emdb_id] = pdb_ids
         else:
              associated_pdbs[emdb_id] = ['-']
         os.remove(temp_xml_file)
         ftp.cwd('../../')
     
    return associated_pdbs

def download_emdb_map(emdb_id):
    ''' To download an emdb map using FTP protocol from the EMDB ID. 
    Returns the downloaded map path 
    
    emdb_id: string. For ex: emdb_id='5778'
    
    returns:
         emdb_map_file. For ex: emdb_map_file='emd_5778.map'
    '''

    EMDB_DIR = "EMD-"+emdb_id 
    EMDB_FILE = "emd_"+emdb_id+'.map.gz'
    
    try:
         command_line_download = "wget ftp://ftp.ebi.ac.uk/pub/databases/emdb/structures/"+EMDB_DIR+'/map/'+EMDB_FILE
         os.system(command_line_download)
    except:
         print("Error downloading EMDB: "+emdb_id)
         return 0
    
    try:
         command_line_extract = "gunzip "+EMDB_FILE
         os.system(command_line_extract)
    except:
         print("Error extracting file EMDB: "+emdb_id)
         return 0
    
    emdb_map_file = "emd_"+emdb_id+".map"
    
    
    if os.path.exists(emdb_map_file):
         EMD_STATUS = "SUCCESS"
         return emdb_map_file,EMD_STATUS
    
    else:
         EMD_STATUS = "FAIL"
         print("Something wrong happened for EMDB "+emdb_id)
         return 0,EMD_STATUS
