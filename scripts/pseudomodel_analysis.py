#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 15 23:54:49 2020

@author: alok
"""

## This script is to analyse Gemmi models saved as PDBs

from pam_headers import *
from emmer.pdb.pdb_tools import *

class Model:
    def __init__(self,points_list):
        self.list = points_list       
        self.unitcell = gemmi.UnitCell(1,1,1,90,90,90)
        self.voxelsize = None
    def calculate_map_values_for_each_point(self,emmap):
        ''' Get map values at each atomic index location.
        
        '''
        
        for point in self.list:
            (x,y,z) = (point.position.x,point.position.y,point.position.z)
            x,y,z = int(round(x)),int(round(y)),int(round(z))
            point.map_value = emmap[z,y,x]
    
    def calculate_nearest_neighbor_dist_for_each_point(self,voxelsize):
        ''' get_neighborhood works only on pixel distance. So use only pixel distance '''
        neighborhood = get_neighborhood(self.list,min_dist_in_pixel=3,only_neighbors=True)
        for i,point in enumerate(self.list):
            point.nearest_neighbor = neighborhood[i][0]*voxelsize
    
    def calculate_relative_acceleration_magnitude(self,emmap,min_dist_in_pixels,g,capmagnitude_map,epsilon,capmagnitude_lj):
        neighborhood = get_neighborhood(self.list,min_dist_in_pixels)
        gz,gy,gx = np.gradient(emmap)
        for i,point in enumerate(self.list):
            lj_neighbors = [self.list[k] for k in neighborhood[i][1]]
            gradient_acceleration,map_value = get_acceleration_from_gradient(gx,gy,gz,emmap, g, point=point, capmagnitude_map=capmagnitude_map)
            
            if len(lj_neighbors)==0:
                lj_potential_acceleration,lj_potential = Vector(np.array([0,0,0])),0
            else:
                lj_potential_acceleration,lj_potential = get_acceleration_from_lj_potential(point, lj_neighbors, epsilon=1, min_dist_in_pixel=min_dist_in_pixels,lj_factor=1.5,capmagnitude_lj=capmagnitude_lj)
            point.gradient_acceleration_magnitude = gradient_acceleration.magnitude()
            point.lj_acceleration_magnitude = lj_potential_acceleration.magnitude()
            try:
                point.relative_acceleration = point.lj_acceleration_magnitude/point.gradient_acceleration_magnitude
            except ZeroDivisionError:
                point.relative_acceleration = 999.99
                print("Zero division error encoutered at: ",+str(point.position.get()))
    
    def get_all_properties(self,emmap,voxelsize,min_dist_in_pixels=1.5,g=10,capmagnitude_map=100,epsilon=1,capmagnitude_lj=100):
        self.calculate_map_values_for_each_point(emmap)
        self.calculate_nearest_neighbor_dist_for_each_point(voxelsize)
        self.calculate_relative_acceleration_magnitude(emmap,min_dist_in_pixels,g,capmagnitude_map,epsilon,capmagnitude_lj)
    
    def extract_pdb_positions(self):
        np_array = np.array([x.pdb_position.get() for x in self.list])
        return np_array
    
    def extract_mrc_positions(self):
        np_array = np.array([x.position.get() for x in self.list])
        return np_array
    
    def get_distance_distribution(self,max_distance):
        np_array = self.extract_pdb_positions()
        sp_tree = spatial.KDTree(np_array)
        
        sparse_distance_matrix= sp_tree.sparse_distance_matrix(sp_tree,max_distance=max_distance)
        return sparse_distance_matrix
    
    def filter_distance(self,min_distance,max_distance):
        distance_matrix = self.get_distance_distribution(max_distance)
        points1 = []
        points2 = []
        for key in distance_matrix.keys():
            if min_distance <= distance_matrix[key] <= max_distance:
                points1.append(key[0])
                points2.append(key[1])
        
        final_points_index = list(set(points1).intersection(set(points2)))
        new_pseudomodel_list = []
        for index in final_points_index:
            point = self.list[index].copy()
            new_pseudomodel_list.append(point)
        new_pseudomodel = AtomList(new_pseudomodel_list)
        return new_pseudomodel
    
    def copy(self):
        new_model = [Atom(x.position.get()) for x in self.list]
        new_model = Model(new_model)
        return new_model
    
    def set_bfactor(self,bfactor):
        for atom in self.list:
            atom.bfactor = bfactor
    
    def combine(self,newmodel):
        set1 = set(self.list)
        set2 = set(newmodel.list)
        combined_list = list(set1.union(set2))
        self.list = []
        self.list = combined_list
       
    def convert_to_gemmi_model(self):
         model = gemmi.Model('pseudo')
         chain_letters = list(string.ascii_uppercase)
         chain_count = 0
         res_count = 0
         atom_count = 0
         model.add_chain(chain_letters[chain_count])
         model = self.add_residue(model,chain_count,res_count)
         for atom in self.list:
             model = self.add_atom(model,chain_count,res_count,atom_count,atom)
             atom_count += 1
             res_count += 1
             model = self.add_residue(model,chain_count,res_count)
                 
             if atom_count % 9999 == 0:
                 chain_count += 1
                 model.add_chain(chain_letters[chain_count])
                 res_count = 0
                 model = self.add_residue(model,chain_count,res_count)
         
         return model
    
    def add_atom(self,model, chain_num, res_num, atom_num, pseudoAtom):
         
         if pseudoAtom.pdb_position.magnitude() == 0:
              
              position = pseudoAtom.position.get()
         else:
              position = pseudoAtom.pdb_position.get()
         atom = gemmi.Atom()
         atom.element = gemmi.Element('O')
         atom.pos = gemmi.Position(position[0],position[1],position[2])
         atom.b_iso = pseudoAtom.bfactor
         atom.occ = pseudoAtom.occ
         
         model[chain_num][res_num].add_atom(atom,atom_num)
         
         return model
        
             
    def add_residue(self,model, chain_num, res_num):
         model[chain_num].add_residue(gemmi.Residue(),res_num)
         model[chain_num][res_num].name = 'HOH'
         model[chain_num][res_num].seqid.num = res_num
     
         return model

    def update_unitcell(self,voxelsize,unitcell=None):
        if unitcell is not None:
            self.unitcell = unitcell
        else:
            voxelsize = self.voxelsize
            num = len(self.list)
            self.unitcell = get_unit_cell_estimate(number_of_atoms=num, vsize=voxel_size)
        self.voxelsize = voxelsize

    def write_pdb(self,output_string,voxelsize,unitcell=None):
          self.update_unitcell(voxelsize,unitcell)
          gemmi_model = self.convert_to_gemmi_model()
          
          structure = gemmi.Structure()
          structure.add_model(gemmi_model)
          structure.cell = self.unitcell
          structure.write_pdb(output_string)
          
    def update_pdb_positions(self,voxelsize=1):
         for atom in self.list:
              atom.pdb_position = atom.position.scale(voxelsize)



def extract_model_from_mask(mask,num_atoms,threshold=1,ignore_these=None):
    from pam_headers import Atom
    
    all_inside_mask = np.asarray(np.where(mask>=threshold)).T.tolist()
    all_inside_set = set([tuple(x) for x in all_inside_mask])
    if ignore_these is not None:
        ignore_set = set([tuple([int(x[0]),int(x[1]),int(x[2])]) for x in ignore_these])
        population = list(all_inside_set.difference(ignore_set))
    else:
        population = list(all_inside_set)
    
    points = [Atom(np.array([x[2]+np.random.uniform(-0.5,0.5),x[1]+np.random.uniform(-0.5,0.5),x[0]+np.random.uniform(-0.5,0.5)])) for x in random.sample(population,num_atoms)]
    model = Model(points)
    return model


    
def get_model_from_gemmi_pdb(pdb_path,emmap_path=None):
    gemmi_model = gemmi.read_pdb(pdb_path)[0]
    
    if emmap_path is not None:
        mrc = mrcfile.open(emmap_path)
        voxelsize = mrc.voxel_size.x
        
        cella = mrc.header.cella
        x = cella.x
        y = cella.y
        z = cella.z
        unitcell = (x,y,z)
    else:
        print(" \n\nWarning! EM-MAP not specified. Setting voxelsize = 1 and unit cell as (1,1,1) \n\n")
        voxelsize = 1
        unitcell = (1,1,1)
    
    points = []
    for chain in gemmi_model:
        for residue in chain:
            for atom in residue:
                position = np.array([atom.pos.x/voxelsize,atom.pos.y/voxelsize,atom.pos.z/voxelsize])
                point = Atom(position)
                point.pdb_position = Vector(position*voxelsize)
                point.bfactor = atom.b_iso
                points.append(point)
    points_list = Model(points)
    points_list.unitcell = unitcell
    
    return points_list

def get_column_array(col_type,points_list):
    ''' 
    Extract point data from the list of points in points_list
    
    '''
    
    column = []
    if col_type == 'relative_acceleration':
        for point in points_list:
            column.append(point.relative_acceleration)
    elif col_type == 'nearest_neighbor_distance':
        for point in points_list:
            column.append(point.nearest_neighbor)
    elif col_type == 'map_value':
        for point in points_list:
            column.append(point.map_value)
    elif col_type == 'bfactor':
        for point in points_list:
            column.append(point.bfactor)
    
    else:
        print("Unknown column heading! ")
    column = np.array(column)
    return column
    

def modify_pdb_column(points,outpdb_name,voxelsize=1,b_factor_type=None,occupancy_type=None):
    '''
    Converts a column, either b factor or occupancy into a desired metric
    
    According to PDB format: b factor and occupancy have max 6 columns, with two decimal places (and one for decimal point)
    Max real value: 999.99
    
    Use voxelsize as 1 if you are extracting from a saved gemmi model
    '''
    b_factor_array = get_column_array(b_factor_type,points.list)  
    occupancy_array = get_column_array(occupancy_type,points.list)
    
    b_factor_array.clip(max=999.99)
    occupancy_array.clip(max=999.99)
    
    gemmi_model = convert_to_gemmi_model(points.list,voxelsize,b_factor_array,occupancy_array)
    write_pdb(gemmi_model,outpdb_name)


                                          
    
    

                
            
        
    
    
        
