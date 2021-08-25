import numpy as np
import string
import gemmi
import matplotlib.pyplot as plt 
import matplotlib.image as image
from time import time
import mrcfile
import random
import math
from scipy.constants import Avogadro
import pandas as pd
from sklearn.neighbors import KDTree
from scipy import spatial
import os,sys
from skimage.color import rgb2gray
from scipy import spatial
from subprocess import call, Popen, PIPE, run
from scipy.signal import correlate
import pandas as pd
from math import pi as pi
from scipy.interpolate import interpn
from emmer.pdb.pdb_tools import *
from datetime import datetime
from pseudomodel_analysis import *

'''
def add_residue(model, chain_num, res_num):
    model[chain_num].add_residue(gemmi.Residue(),res_num)
    model[chain_num][res_num].name = 'HOH'
    model[chain_num][res_num].seqid.num = res_num

    return model

def add_atom(model, chain_num, res_num, atom_num,position,voxelsize,b_factor=20,occupancy=1):
    if len(position) == 2:
        position = np.array([position[0],position[1],0])
    atom = gemmi.Atom()
    atom.element = gemmi.Element('O')
    atom.pos = gemmi.Position(position[0]*voxelsize,position[1]*voxelsize,position[2]*voxelsize)
    atom.b_iso = b_factor
    atom.occ = occupancy
    
    model[chain_num][res_num].add_atom(atom,atom_num)

ci
def write_pdb(pseudomodel,output_string):
    structure = gemmi.Structure()
    structure.add_model(pseudomodel)
    structure.write_pdb(output_string)

def convert_to_gemmi_model(points,voxelsize,b_factor_array=None,occupancy_array=None):
    np_points = np.array([x.position.get() for x in points])
    model = gemmi.Model('pseudo')
    chain_letters = list(string.ascii_uppercase)
    chain_count = 0
    res_count = 0
    atom_count = 0
    model.add_chain(chain_letters[chain_count])
    model = add_residue(model,chain_count,res_count)
    if b_factor_array is None:
        b_factor_array = np.ones(len(points))*20
    if occupancy_array is None:
        occupancy_array = np.ones(len(points))
    for ix,point in enumerate(np_points):
        model = add_atom(model,chain_count,res_count,atom_count,point,voxelsize,b_factor_array[ix],occupancy_array[ix])
        atom_count += 1
        res_count += 1
        model = add_residue(model,chain_count,res_count)
            
        if atom_count % 9999 == 0:
            chain_count += 1
            model.add_chain(chain_letters[chain_count])
            res_count = 0
            model = add_residue(model,chain_count,res_count)
    
    return model

'''


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

def value_at_point(g,x,y,z):
    return g[z,y,x]

def get_gradient(g,point):
    [xi,yi,zi] = [int(round(point.position.x)),int(round(point.position.y)),int(round(point.position.z))]
    x = np.arange(xi-1,xi+2,1)
    y = np.arange(yi-1,yi+2,1)
    z = np.arange(zi-1,zi+2,1)
    points = (z,y,x)
    values = g[zi-1:zi+2,yi-1:yi+2,xi-1:xi+2]
  #  print(points[0].shape)
  #  print(values.shape)
    (xr,yr,zr) = point.position.get()
    interpolation_point = np.array([zr,yr,xr])
    interpolated_value = interpn(points,values,interpolation_point,method='linear')
    
    return interpolated_value[0]    

        
def get_acceleration_from_gradient(gx,gy,gz,emmap,g,point,capmagnitude_map):
    [x,y,z] = [int(round(point.position.x)),int(round(point.position.y)),int(round(point.position.z))]

    
    theta_x = gx[z,y,x] 
    theta_y = gy[z,y,x] 
    theta_z = gz[z,y,x] 
    '''
    
    theta_x = get_gradient(gx,point)
    theta_y = get_gradient(gy,point)
    theta_z = get_gradient(gz,point)
    
    '''
    
    acceleration_x = g * theta_x
    acceleration_y = g * theta_y
    acceleration_z = g * theta_z

    
    acceleration = Vector(np.array([acceleration_x,acceleration_y,acceleration_z]))

    
    return acceleration.cap_magnitude(capmagnitude_map),emmap[z,y,x]


def get_acceleration_from_lj_potential(targetpoint,lj_neighbors,epsilon,min_dist_in_pixel,lj_factor,capmagnitude_lj):
    #print(lj_factor)
    tic = time()

    lj_neighbors_points = [x.position.get() for x in lj_neighbors]
    distance_vector = targetpoint.position.get() - lj_neighbors_points
    r = np.sqrt(np.einsum('ij->i',distance_vector**2))
    unit_diff_vector = (distance_vector.transpose() / r).transpose()

    
    eps = epsilon
    rm = min_dist_in_pixel*lj_factor
    rm1 = 1.23
    rm2 = 1.33
    rm3 = 1.45
    rm4 = 1.52
  
    v_lj = eps * ((rm1/r)**12 - 2*(rm1/r)**6)
    
    f_r = np.array((12 * eps * rm**6 * (r**6 - rm**6))/r**13)
    #f_r = np.array((288 * eps * rm1**48 * (r**48 - rm1**48))/r**97) + np.array((288 * eps * rm2**48 * (r**48 - rm2**48))/r**97) + np.array((288 * eps * rm3**48 * (r**48 - rm3**48))/r**97) + np.array((288 * eps * rm4**48 * (r**48 - rm4**48))/r**97)
    '''
    f_r = []
    for i in range(len(r)):
        if 0 <= r[i] < 1.3:
            f_r.append(288 * eps * rm1**48 * (r[i]**48 - rm1**48)/r[i]**97)
        elif 1.3 <= r[i] < 1.4:
            f_r.append(288 * eps * rm2**48 * (r[i]**48 - rm2**48)/r[i]**97)
        elif 1.4 <= r[i] < 1.5:
            f_r.append(288 * eps * rm3**48 * (r[i]**48 - rm3**48)/r[i]**97)
        else:
            f_r.append(288 * eps * rm4**48 * (r[i]**48 - rm4**48)/r[i]**97)
    
    f_r = np.array(f_r)
    '''
    
    f_r_vector = np.array([np.array(f_r[k])*np.array(unit_diff_vector[k]) for k in range(len(lj_neighbors))])
    
  
    fx = -f_r_vector[:,0]
    fy = -f_r_vector[:,1]
    fz = -f_r_vector[:,2]

    
    ax = fx.sum() / targetpoint.mass
    ay = fy.sum() / targetpoint.mass
    az = fz.sum() / targetpoint.mass
    acc = Vector(np.array([ax,ay,az]))


    return acc.cap_magnitude(capmagnitude_lj),v_lj.sum()

def get_neighborhood(points,min_dist_in_pixel,fromArray=False,only_neighbors=False):
    '''
    input: points is a list of point objects. If the list is already a numpy array then the variable fromArray must be True
    rerturn a dictionary of neighborhoods. If only_neighbors is true, then only distance to neighbor is sent. Else, distance to neighbor and the indices of all nearest neighbors (distance of min_dist * 3 ) is sent
    '''
    
    tic = time()
    if fromArray==False:
        np_points = np.array([list(x.position.get()) for x in points])
    else:
        np_points = points
    neighborhood = {}
    tree = KDTree(np_points)
    if only_neighbors==False:
        for i in range(len(points)):
            ind = tree.query_radius(np_points[i:i+1],r=min_dist_in_pixel*3)[0]
            d,ix = tree.query(np_points[i:i+1],k=2)
            ind = np.delete(ind,np.where(ind==i))
            neighborhood[i]=[d[0][1],ind]
        return neighborhood
    else:
        for i in range(len(points)):
            d,ix = tree.query(np_points[i:i+1],k=2)
            neighborhood[i]=[d[0][1]]
        return neighborhood
        


def average_map_value(points):
    sum_map_value = 0
    for point in points:
        sum_map_value += point.map_value
    average_mapvalue = sum_map_value/len(points)
    return average_mapvalue

def compute_real_space_correlation(map1,map2):
    '''
    map1 and map2 are both numpy arrays with same of shape (H,B,W)
    output: RSCC which compares the cross correlation (normalised) between the two maps
    
    '''
    (map1_mean,map1_std) = (map1.mean(),map1.std())
    (map2_mean,map2_std) = (map2.mean(),map2.std())
    
    n = map1.size
    
    RSCC = (((map1-map1_mean)*(map2-map2_mean))/(map1_std*map2_std)).sum() * (1/n)
    
    return RSCC

def main_solver3D(emmap,gx,gy,gz,model_initial,g,friction,min_dist_in_angst,voxelsize,
                  dt=0.05,capmagnitude_lj=400,epsilon=1,scale_lj=1,lj_factor=1,capmagnitude_map=100,scale_map=1,total_iterations=50, 
                  path_for_gemmi_models=None,emmap_path=None,mask_path=None,returnPointsOnly=True,verbose=False,
                  integration='verlet',myoutput=None):
    '''
    Function to solve pseudoatomic model using gradient descent approach. 
    
    emmap : numpy.ndarray
        Numpy array containing the 3D volume of the map
    gx,gy,gz : numpy.ndarray
        Gradients obtained using numpy.gradient() method to get gradient information in x,y and z
    model_initial : pseudomodel_analysis.Model()
        Is a custom built class which has the coordinate information of all atoms. Also has several useful custom functions 
    g : float
        Gradient scaling parameter to scale the "accelerations" uniformly across the model
    friction : float
        friction coefficient to converge the model
    min_dist_in_angst : float
        Minimum distance between two atoms in the pseudo-atomic model, constrained by the bond lengths
    voxelsize : float
        Voxelsize of the emmap
    
    -- special note for the following parameters --
    capmagnitude_lj, capmagnitude_map : float
        These values truncate the maximum acceleration felt by an atom during each iteration so that the analysis becomes bounded
    
        
    '''
    
    peak_bond_length_list = []
    map_values = []
    pseudomodel = model_initial.copy()
    gradient_magnitude = np.sqrt(gx**2+gy**2+gz**2)
    if verbose:
        solver_properties = 'Solver started with the following properties: \n'+'\n Number of atoms = '+str(len(pseudomodel.list))+'\n Map potential: \n'+'\n g = '+str(g)+'\n Max gradient magnitude  = '+str(gradient_magnitude.max())+'\n Map value range  = '+str((emmap.min(),emmap.max()))+'\n Cap magnitude at  = '+str(capmagnitude_map)+'\n LJ Potential: \n'+'\n Equilibrium distance = '+str(min_dist_in_angst)+'\n Voxelsize, in A = '+str(voxelsize)+'\n LJ Factor = '+str(lj_factor)+'\n Epsilon = '+str(epsilon)+'\n Cap magnitude at  = '+str(capmagnitude_lj)+'\n Friction: \n'+ '\n Friction Coefficient = '+str(friction)+'\n Solver properties: \n'+'\n Total Iterations = '+str(total_iterations)+'\n Time step = '+str(dt)
              
        print(solver_properties)
        if myoutput is not None:
            myoutput.write(solver_properties)
    if verbose:    
        print('# \t|\t Peak bond length \t | \t Minimum bond length \t | \t Average map value')        
    for iter in range(total_iterations):
        
        neighborhood = get_neighborhood(pseudomodel.list,min_dist_in_angst/voxelsize)
        all_bond_lengths = np.array([d[0]*voxelsize for d in neighborhood.values()])
        bond_length_histogram = np.histogram(all_bond_lengths,bins=200)
        peak_bond_length = bond_length_histogram[1][bond_length_histogram[0].argmax()]
        peak_bond_length_list.append(peak_bond_length)
    
        point_id = 0
        for atom in pseudomodel.list:
            lj_neighbors = [pseudomodel.list[k] for k in neighborhood[point_id][1]]
            
            gradient_acceleration,map_value = get_acceleration_from_gradient(gx,gy,gz,emmap, g, point=atom, capmagnitude_map=capmagnitude_map)
            if len(lj_neighbors)==0:
                lj_potential_acceleration,lj_potential = Vector(np.array([0,0,0])),0
            else:
                lj_potential_acceleration,lj_potential = get_acceleration_from_lj_potential(atom, lj_neighbors, epsilon=1, min_dist_in_pixel=min_dist_in_angst/voxelsize,lj_factor=lj_factor,capmagnitude_lj=capmagnitude_lj)
            # TO BE CONTINUES
            gradient_acceleration,lj_potential_acceleration = gradient_acceleration.scale(scale_map),lj_potential_acceleration.scale(scale_lj)
            acceleration = add_Vector(gradient_acceleration,lj_potential_acceleration)
            # add friction 
            atom.acceleration = add_Vector(acceleration, atom.velocity.scale(-friction))
            atom.map_value = map_value
            point_id += 1
        
        map_values.append(average_map_value(pseudomodel.list))
        
        if integration == 'euler':
            for atom in pseudomodel.list:
                atom.velocity_from_acceleration(dt)        
                atom.position_from_velocity(dt)
                atom.update_history()
        
        elif integration == 'verlet':
            ''' 
            For the first iteration, use Euler integration since we have no information about -1'th time step
            ''' 
            if iter == 0: 
                for atom in pseudomodel.list:
                    atom.velocity_from_acceleration(dt)        
                    atom.position_from_velocity(dt)
                    atom.update_history()
            else:
                for atom in pseudomodel.list:
                    atom.verlet_integration(dt)
                    atom.update_history()
        
            
        
        
        if path_for_gemmi_models is not None:
            gemmi_model = convert_to_gemmi_model(pseudomodel.list,voxelsize);
            
            emmap_shape = emmap.shape
            unitcell = gemmi.UnitCell(emmap_shape[0]*voxelsize,emmap_shape[1]*voxelsize,emmap_shape[2]*voxelsize)
            
            gemmi_structure = gemmi.Structure()
            gemmi_structure.add_model(gemmi_model)
            gemmi_structure.cell = unitcell
            
            map_iteration = pdb_to_map(pdb_structure=gemmi_structure,vsize=voxelsize)
            size = map_iteration.size
            cc = correlate(normalize(map_iteration), normalize(emmap),mode='same',method='fft') / size
            cross_correlation.append(cc.max())
    
    
            if verbose: 
                print(str(iter)+": #peak_bond_length = "+str(peak_bond_length_list[iter])+": #map_value = "+str(map_values[iter])+": cc_max = "+str(cc.max()))    
    
        else:
            if verbose:
                print(str(iter)+
                      "\t | \t "+str(peak_bond_length_list[iter])+
                      "\t | \t "+str(all_bond_lengths.min())+
                      "\t | \t "+str(map_values[iter]))    
            
    
    if returnPointsOnly == True:
        return pseudomodel    
    else:
        if path_for_gemmi_models is not None:
            return pseudomodel, peak_bond_length_list, map_values, cross_correlation
        else:
            return pseudomodel, peak_bond_length_list, map_values


def find_and_kick(points_array,kicklist,kick):
    '''
    Function to return a disturbed point cloud, given an input point cloud

    Parameters
    ----------
    points_array : numpy.ndarray
        A numpy array, where each element in the array has a coordinate information. 
    kicklist : list
        Index of atoms which need to be "kicked" by a random value between +kick and -kick
    kick : int
        Magnitude of "kick" to an atom in a given direction

    Returns
    -------
    points_array : numpy.ndarray
        A numpy array, where each element in the array has a coordinate information of disturbed point clouds

    '''
    N = len(kicklist)
    points_array = np.array(points_array,dtype=float)
    for i in range(N):
        points_array[kicklist[i]]+=[random.uniform(-kick,kick),random.uniform(-kick,kick),random.uniform(-kick,kick)]
    return points_array    

        
def main_solver_kick(model_initial, min_dist_in_angst, voxelsize, total_iterations=99,returnPointsOnly=True,verbose=False):
    '''
    Solver to iteratively morph a point cloud so that it satisfies a minimum distance criterion for any pair of points

    Parameters
    ----------
    model_initial : pseudomodel_analysis.Model()
        Is a custom built class which has the coordinate information of all atoms before satisfying minimum distance criteria
    min_dist_in_angst : float
        Minimum distance between two atoms in the pseudo-atomic model, constrained by the bond lengths
    voxelsize : float
        Voxelsize of the emmap
    total_iterations : int, optional
        Number of iterations to run the solver. The default is 99.
    returnPointsOnly : bool, optional
        If true, returns only the model. If false, it returns other analysis parameters. The default is True.
    
    Returns
    -------
    pseudomodel : pseudomodel_analysis.Model()
        Is a custom built class which has the coordinate information of all atoms after satisfying minimum distance criteria 

    '''
    points_array = np.array([x.position.get() for x in model_initial.list])
    number_of_contacts = []
    if verbose:
        print(' Solver started with the following properties: \n'+
              '\n Number of atoms = '+str(len(points_array))+
              '\n Equilibrium distance = '+str(min_dist_in_angst)+
              '\n Voxelsize, in A = '+str(voxelsize)+
              '\n Total Iterations = '+str(total_iterations))
              
    for i in range(total_iterations):
        neighbors = get_neighborhood(points_array,min_dist_in_angst,fromArray=True)
        kicklist = [x for x in neighbors.keys() if neighbors[x][0] <= min_dist_in_angst/voxelsize]
        points_array = find_and_kick(points_array,kicklist,kick=1)
        number_of_contacts.append(len(kicklist))
        if verbose: 
            print("Iteration number =  "+str(i)+": # Atoms less than eq. dist = "+str(len(kicklist)))
        
        if sum(number_of_contacts[-3:]) == 0:
            break
        
    pseudomodel = Model([Atom(x) for x in points_array])
    pseudomodel.voxelsize = voxelsize
    pseudomodel.update_pdb_positions(voxelsize)
    if returnPointsOnly:
        return pseudomodel
    else:
        return pseudomodel, number_of_contacts

def compute_radial_profile(vol, center=[0,0,0], return_indices=False):
    dim = vol.shape
    m = np.mod(vol.shape,2)
    # make compliant with both fftn and rfftn
    if center is None:
        ps = np.abs(np.fft.fftshift((np.fft.fftn(vol))))
        z, y, x = np.indices(ps.shape)
        center = tuple((a - 1) / 2.0 for a in ps.shape[::-1])
        radii = np.sqrt((x - center[0])**2 + (y - center[1])**2 + (z - center[2])**2)
        radii = radii.astype(np.int)
    else:
        ps = np.abs( np.fft.rfftn(vol) )
        if not return_indices:
            x, y, z = np.indices(ps.shape)
            radii = np.sqrt(x**2 + y**2 + z**2)
            radii = radii.astype(np.int)
        else:
            [x, y, z] = np.mgrid[-dim[0]//2+m[0]:(dim[0]-1)//2+1, -dim[1]//2+m[1]:(dim[1]-1)//2+1, 0:dim[2]//2+1]
            x = np.fft.ifftshift(x)
            y = np.fft.ifftshift(y)
            radii = np.sqrt(x**2 + y**2 + z**2)
            radii = radii.astype(np.int)
    radial_profile = np.bincount(radii.ravel(), ps.ravel()) / np.bincount(radii.ravel())
    # exclude corner frequencies
    radial_profile = radial_profile[0:int(ps.shape[0]/2)]
    if not return_indices:
        return radial_profile
    else:
        return radial_profile, radii
    
def compute_and_plot_radial_profiles_from_gemmi_pdb(paths_to_pdb,path_to_refmap_script,emmap_path,mask_path,apix):
    colors = ['r','g','b','k','w']
    radial_profiles = []
    maps = []
    for i,path in enumerate(paths_to_pdb):
        command_line = "phenix.python "+path_to_refmap_script+" -mc "+path+" -em "+emmap_path+" -ma "+mask_path
        os.system(command_line)
        modelmap_path = path[:-4]+'_4locscale.mrc'
        model_map = mrcfile.open(modelmap_path).data
        radial_profile = compute_radial_profile(model_map)
        freq = np.linspace(1./(float(apix)*radial_profile.shape[0]), 1./(float(apix)*2), radial_profile.shape[0],endpoint=True)
        #plt.subplot(1,2,1)
        #plt.plot(freq,radial_profile/radial_profile.max(),colors[i])
        #plt.subplot(1,2,2)2156
        plt.plot(freq[30:],radial_profile[30:]/radial_profile.max(),colors[i])
        radial_profiles.append(radial_profile)
        maps.append(model_map)
    return maps,radial_profiles


def compute_radial_profile_from_mrcs(mrc_paths,keys=None,slice_index=30):
    if keys is None:
        keys = mrc_paths
    
    mrcs = []
    for mrc in mrc_paths:
        mrcs.append(mrcfile.open(mrc))
    
    k = 0
    emmaps = {}
    radial_profiles = {}
    freq={}
    
    for mrc in mrcs:
        emmaps[keys[k]] = mrc.data
        k += 1
    for key in keys:
        radial_profiles[key] = compute_radial_profile(emmaps[key])
        radial_profiles[key] /= radial_profiles[key].max()
    k = 0
    for key in keys:
        shapes = radial_profiles[key].shape[0]
        apix = mrcs[k].voxel_size.x
        freq[key] = np.linspace(1./(float(apix)*shapes), 1./(float(apix)*2), shapes,endpoint=True)
        k += 1 
    
    for key in keys:
        plt.plot(freq[key][slice_index:],radial_profiles[key][slice_index:])
    plt.legend(keys)
    
    for key in keys:
        radial_profiles[key] = tuple([freq,radial_profiles[key]])
    return radial_profiles, emmaps
    
        

def add_cryst1_line(pdb_path,unitcell=None,emmap_path=None,new_pdb_path=None):
    '''
    pdb_path -> Address of .pdb path
         
    Some PDB files developed for cryoEM maps do not have proper cryst1 record. Two options to modify:

    1. From an input tuple, or array. In this case, unitcell is a python tuple, which has unit cell dimensions in angstorm
    Ex: unitcell = (x,y,z)
    2. From a mrcfile. In this case, point to an associated EM map and the unit cell dimensions are taken from that
    emmap_path -> Address of associated .mrc file
    
    If you like to save the pdb file with a different name, or address then change the 'new_pdb_path' 
    
    '''
    if emmap_path is not None:
        mrc = mrcfile.open(emmap_path)
        cella = mrc.header.cella
        x = cella.x
        y = cella.y
        z = cella.z
    elif unitcell is not None:
        x = unitcell[0]
        y = unitcell[1]
        z = unitcell[2]
    else:
        print("Please give either unit cell dimensions (in Ang) or point to an associated mrc file!")
        return
    
    unitcell = gemmi.UnitCell(x,y,z,90,90,90)
    
    gemmi_structure = gemmi.read_pdb(pdb_path)
    gemmi_structure.cell = unitcell
    if new_pdb_path is None:
        gemmi_structure.write_pdb(pdb_path)
    else:
        gemmi_structure.write_pdb(new_pdb_path)
    