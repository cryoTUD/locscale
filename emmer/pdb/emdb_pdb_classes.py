#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 24 14:52:38 2021

@author: alok
"""
class date:
     def __init__(self,date_string,date_format='yyyy-mm-dd'):
          datestring = date_string.split('-')
          if date_format == 'yyyy-mm-dd':
               self.year = int(datestring[0])
               self.month = int(datestring[1])
               self.date = int(datestring[2])              
               
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
               print("Default values set. Now updating the properties by downloading and extracting the XML file")
               self.parse_properties()
          
          
     def download(self):
               
          ''' To download an emdb map using FTP protocol from the EMDB ID. 
          Returns the downloaded map path 
          
          emdb_id: string. For ex: emdb_id='5778'
          
          returns:
               emdb_map_file. For ex: emdb_map_file='emd_5778.map'
          '''
          import os
          EMDB_DIR = "EMD-"+self.id
          EMDB_FILE = "emd_"+self.id+'.map.gz'
          emdb_id = self.id
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
               print("Map downloaded! You can see the file here: "+emdb_map_file)
               return emdb_map_file
          
          else:
               
               print("Something wrong happened for EMDB "+emdb_id)
               return 0
          
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
          print("Now downloading the XML file")
          header_directory = 'EMD-'+emdb_id+'/header'
          xml_file = 'emd-'+emdb_id+'.xml'
          ftp.cwd(header_directory)
          ftp.retrbinary("RETR "+xml_file,open(filepath,'wb').write)
          if os.path.exists(filepath):
               print("XML downloaded! You can find the XML file here: "+filepath)
               return filepath
          else:
               print("Coudl not download the XML file :( ")
               return 0
     
     def parse_properties(self,header_file=None,delete_xml=True):
          from xml.dom import minidom
          import os
          #import xml.etree.ElementTree as ET
          
          if header_file is None:
               print("You did not pass a header file. Now downloading from the server!")
               header_file = self.get_header_file()
          
          xmldoc = minidom.parse(header_file)
          #tree = ET(header_file)
                    
          self.header_info_xml = xmldoc
          self.parse_deposition()
          self.parse_processing()
          if delete_xml:
               print("Deleting XML file. Set delete_xml=False to disable this!")
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
               print("Problem with finding fitted PDBS")
               
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
     
     

          