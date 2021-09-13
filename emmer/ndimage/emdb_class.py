#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 24 14:52:38 2021

@author: alok
"""

import os
import gzip
import shutil

class date:
    def __init__(self,date_string,date_format='yyyy-mm-dd'):
         datestring = date_string.split('-')
         if date_format == 'yyyy-mm-dd':
              self.year = int(datestring[0])
              self.month = int(datestring[1])
              self.date = int(datestring[2])  
    def show(self):
        print(str(self.year) + "-" + str(self.month) + "-" + str(self.date))
    def get(self):
        return str(self.year) + "-" + str(self.month) + "-" + str(self.date)
               
class EMDB:
     '''
     * Create a class called emdb and define the properties
     * Degine methods to retrive the properties of an emdb object
     '''
     def __init__(self,emdbid,parse_properties=True):
          # emdbid is a string. For example: emdbid = '5778'
          if isinstance(emdbid,str):
               self.id = emdbid
          else:
               self.id = str(emdbid)
          self.resolution = 0
          #self.pixelsize = np.array([0,0,0])
          self.fitted_pdbs = ['']
          self.deposited_date = 0
          self.header_info_xml = 0
          self.method = 'unknown'
          if parse_properties:
               #print("Default values set. Now updating the properties" +
               #      "by downloading and extracting the XML file")
               self.parse_properties()
          
          
     def download(self, dataset_dir="EMDB_maps", num_of_nested_dirs=1):
               
          ''' To download an emdb map using FTP protocol from the EMDB ID. 
          Returns the downloaded map path 
          
          emdb_id: string. For ex: emdb_id='5778'
          num_of_nested_dirs: int. For ex: 1, 2, 3 - dependent on dataset_dir structure
          
          returns:
                 emdb_map_file. For ex: emdb_map_file='emd_5778.map'
          '''
          
          emdb_id = self.id
          EMDB_DIR = "EMD-" + emdb_id
          EMDB_FILE = "emd_" + emdb_id + '.map.gz'
          
          try:
               os.chdir(dataset_dir)
          except FileNotFoundError:
               os.mkdir(dataset_dir)
               os.chdir(dataset_dir)
          
          command_line_download = "wget --quiet ftp://ftp.ebi.ac.uk/pub/"\
               + "databases/emdb/structures/" + EMDB_DIR \
               + '/map/' + EMDB_FILE
          os.system(command_line_download)
          
          with gzip.open("emd_" + emdb_id + ".map.gz", 'rb') as f_in:
               with open("emd_" + emdb_id + ".map", 'wb') as f_out:
                    shutil.copyfileobj(f_in, f_out)
          for file in os.listdir():
               if file.startswith("emd_" + emdb_id + ".map.gz"):
                    os.remove(file)
          for i in range(num_of_nested_dirs):
               os.chdir("..")
          print("- Successfully downloaded EMDBmap {} from database".format(emdb_id))
          
     def get_header_file(self,filepath=None):
          if filepath is None:
               filepath = self.id + '.xml'

          import ftplib
          import os
                    
          emdb_id = self.id
          
          ftp_host = "ftp.ebi.ac.uk"
          ftp = ftplib.FTP(ftp_host,'anonymous')
          directory = "pub/databases/emdb/structures/"
          ftp.cwd(directory)
          # print("Now downloading the XML file")
          header_directory = 'EMD-'+emdb_id+'/header'
          xml_file = 'emd-'+emdb_id+'.xml'
          ftp.cwd(header_directory)
          ftp.retrbinary("RETR "+xml_file,open(filepath,'wb').write)
          if os.path.exists(filepath):
               # print("XML downloaded! You can find the XML file here: "+filepath)
               return filepath
          else:
               print("Coudl not download the XML file :( ")
               return 0
     
     def parse_properties(self,header_file=None,delete_xml=True):
          from xml.dom import minidom
          import os
          #import xml.etree.ElementTree as ET
          
          if header_file is None:
               # print("You did not pass a header file. Now downloading from the server!")
               header_file = self.get_header_file()
          
          xmldoc = minidom.parse(header_file)
          #tree = ET(header_file)
                    
          self.header_info_xml = xmldoc
          self.parse_deposition()
          self.parse_processing()
          if delete_xml:
               # print("Deleting XML file. Set delete_xml=False to disable this!")
               os.remove(header_file)
          

     def parse_deposition(self):
          xmlfile = self.header_info_xml
          try:
               deposition = xmlfile.getElementsByTagName('deposition')[0]
               deposition_date = deposition.getElementsByTagName('depositionDate')[0].childNodes[0].nodeValue
               self.deposited_date = date(deposition_date)
          except:
               print("Problem with deposited date")
          
          try:
               fitted_pdbs = deposition.getElementsByTagName('fittedPDBEntryIdList')[0]
               fitted_pdb_list = []
               for pdb_entry in fitted_pdbs.getElementsByTagName('fittedPDBEntryId'):
                    pdb_id = pdb_entry.childNodes[0].nodeValue
                    fitted_pdb_list.append(pdb_id)
               
               self.fitted_pdbs = fitted_pdb_list
          except:
               pass
               #print("Problem with finding fitted PDBS")
               
     def parse_processing(self):
          xmlfile = self.header_info_xml
          try:
               processing = xmlfile.getElementsByTagName('processing')[0]
               self.method = processing.getElementsByTagName('method')[0].childNodes[0].nodeValue
          except:
               print("Problem finding method")
          
          try:
               reconstruction = processing.getElementsByTagName('reconstruction')[0]
               self.resolution = reconstruction.getElementsByTagName('resolutionByAuthor')[0].childNodes[0].nodeValue
          except:
               print("Problem finding resolution")
               
def get_all_emdb_id():
     import ftplib
     import os
     
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
     
     import ftplib
     import os
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