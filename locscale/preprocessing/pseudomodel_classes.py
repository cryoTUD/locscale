#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 15 23:54:49 2020

@author: alok
"""

## This script is to analyse Gemmi models saved as PDBs
import numpy as np
import math
import gemmi
import mrcfile
#from pam_headers import *
#from emmer.pdb.pdb_tools import *


class Vector:
    def __init__(self,input_array):
        self.x = input_array[0]
        self.y = input_array[1]
        self.z = input_array[2]
    def get(self):
        return np.array([self.x,self.y,self.z])
    def magnitude(self):
        return math.sqrt(self.x**2+self.y**2+self.z**2)
    def cap_magnitude(self,cap):
        mag = self.magnitude()
        if mag > cap:
            factor= cap/mag
            return self.scale(factor)
        else:
            return self
    
    def scale(self,scale):
        return Vector(scale*self.get())

def add_Vector(vector_a,vector_b):
        return Vector(vector_a.get() + vector_b.get())

d_type = [('pos',tuple),('vel',tuple),('acc',tuple)]

class Atom:
    def __init__(self,init_pos):
        self.id = 0
        self.position = Vector(init_pos)
        self.pdb_position = Vector(np.array([0,0,0]))
        self.velocity = Vector(np.array([0,0,0]))    
        self.acceleration = Vector(np.array([0,0,0]))
        self.mass = 1 # Mass factor - not in real units! 
        self.nearest_neighbor = math.inf
        self.gradient_acceleration_magnitude = 0
        self.lj_acceleration_magnitude = 0
        self.relative_acceleration = 0
        self.map_value = 0
        self.bfactor = 20
        self.occ = 1
        self.position_history = [self.position.get()]
        self.velocity_history = [self.velocity.get()]
        self.acceleration_history = [self.acceleration.get()]
        self.map_value_history = [self.map_value]
        
    def get_distance_vector(self,target):
        distance_vector = Vector(np.array(self.position.get() - target.position.get()))
        return distance_vector
    
    def angle_wrt_horizontal(self,target):
        return math.atan2(target.position.y - self.position.y, target.position.x - self.position.x)
    
    def velocity_from_acceleration(self,dt):
        vx = self.velocity.x + self.acceleration.x*dt
        vy = self.velocity.y + self.acceleration.y*dt
        vz = self.velocity.z + self.acceleration.z*dt
       # print('velocity: '+str(tuple([vx,vy])))
        self.velocity = Vector(np.array([vx,vy,vz]))
        
    def position_from_velocity(self,dt):
        x = self.position.x + self.velocity.x*dt
        y = self.position.y + self.velocity.y*dt
        z = self.position.z + self.velocity.z*dt
        self.position = Vector(np.array([x,y,z]))
    def verlet_integration(self,dt):
        ## Update positions
        
        r_now = self.position_history[-1]
        r_prev = self.position_history[-2]
        a_now = self.acceleration.get()
        
        r_next = 2 * r_now - r_prev + a_now * dt**2
        
        self.position = Vector(r_next)
        
        # Update velocities 
        
        v_next = (r_next - r_prev) / (2*dt)
        
        self.velocity = Vector(v_next)
        
    def update_history(self):
        self.position_history.append(self.position.get())
        self.velocity_history.append(self.velocity.get())
        self.acceleration_history.append(self.acceleration.get())
        self.map_value_history.append(self.map_value)
    
    def copy(self):
        position = self.position.get()
        newPoint = Atom(position)
        newPoint.pdb_position = self.pdb_position
        return newPoint

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
        from locscale.preprocessing.pseudomodel_solvers import get_neighborhood
        
        neighborhood = get_neighborhood(self.list,min_dist_in_pixel=3,only_neighbors=True)
        for i,point in enumerate(self.list):
            point.nearest_neighbor = neighborhood[i][0]*voxelsize
    
    def calculate_relative_acceleration_magnitude(self,emmap,min_dist_in_pixels,g,capmagnitude_map,epsilon,capmagnitude_lj):
        
        from locscale.preprocessing.pseudomodel_solvers import get_neighborhood, get_acceleration_from_gradient, get_acceleration_from_lj_potential
        
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
        
        from scipy import spatial
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
        new_pseudomodel = Model(new_pseudomodel_list)
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
        import string
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
         atom.name = "O"
         
         model[chain_num][res_num].add_atom(atom,atom_num)
         
         return model
        
             
    def add_residue(self,model, chain_num, res_num):
         model[chain_num].add_residue(gemmi.Residue(),res_num)
         model[chain_num][res_num].name = 'HOH'
         model[chain_num][res_num].seqid.num = res_num
     
         return model

    def update_unitcell(self,voxelsize,unitcell=None):
        from locscale.include.emmer.pdb.pdb_tools import get_unit_cell_estimate
        
        if unitcell is not None:
            self.unitcell = unitcell
        else:
            voxelsize = self.voxelsize
            num = len(self.list)
            self.unitcell = get_unit_cell_estimate(number_of_atoms=num, vsize=voxelsize)
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
    from locscale.preprocessing.pseudomodel_classes import Atom
    import random
    x1,x2,x3 = mask.shape
    buffer = 2 ## To ensure no atoms near edge get picked 
    ones_array = np.ones((x1-2*buffer, x2-2*buffer, x3-2*buffer))
    padded = np.pad(ones_array, buffer)
    
    edge_cropped_mask = mask * padded
    
    
    all_inside_mask = np.asarray(np.where(edge_cropped_mask>=threshold)).T.tolist()
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
    



                                          
    
    

                
            
        
    
    
        
