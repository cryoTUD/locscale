import numpy as np
       
def get_acceleration_from_gradient(gx,gy,gz,emmap,g,point,capmagnitude_map):
    from locscale.preprocessing.pseudomodel_classes import Vector
    
    [x,y,z] = [int(round(point.position.x)),int(round(point.position.y)),int(round(point.position.z))]
    
    try:
    
        theta_x = gx[z,y,x] 
        theta_y = gy[z,y,x] 
        theta_z = gz[z,y,x] 
    except Exception as e:
        print((x,y,z))
        print(e)
        raise

    
    acceleration_x = g * theta_x
    acceleration_y = g * theta_y
    acceleration_z = g * theta_z

    
    acceleration = Vector(np.array([acceleration_x,acceleration_y,acceleration_z]))

    
    return acceleration.cap_magnitude(capmagnitude_map),emmap[z,y,x]


def get_acceleration_from_lj_potential(targetpoint,lj_neighbors,epsilon,min_dist_in_pixel,lj_factor,capmagnitude_lj):
    from locscale.preprocessing.pseudomodel_classes import Vector

    lj_neighbors_points = [x.position.get() for x in lj_neighbors]
    distance_vector = targetpoint.position.get() - lj_neighbors_points
    r = np.sqrt(np.einsum('ij->i',distance_vector**2))
    unit_diff_vector = (distance_vector.transpose() / r).transpose()

    
    eps = epsilon
    rm = min_dist_in_pixel*lj_factor

    v_lj = eps * ((rm/r)**12 - 2*(rm/r)**6)
    
    f_r = np.array((12 * eps * rm**6 * (r**6 - rm**6))/r**13)
        
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
    from sklearn.neighbors import KDTree

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
    map_val = []
    for point in points:
        map_val.append(point.map_value)
    map_val = np.array(map_val)
    
    average_mapvalue = round(map_val.mean(),3)
    sd_mapvalue = round(map_val.std(),3)
    return (average_mapvalue,sd_mapvalue)

def acceleration_contribution(pseudomodel):
    gradient_list= []
    lj_list= []
    for atom in pseudomodel.list:
        gradient_list.append(atom.gradient_acceleration_magnitude)
        lj_list.append(atom.lj_acceleration_magnitude)
    
    return gradient_list, lj_list
    


def main_solver3D(emmap,gx,gy,gz,model_initial,g,friction,min_dist_in_angst,voxelsize,
                  dt=0.05,capmagnitude_lj=400,epsilon=1,scale_lj=1,lj_factor=1,capmagnitude_map=100,scale_map=1,total_iterations=50, 
                  compute_map=False,emmap_path=None,mask_path=None,returnPointsOnly=True,verbose=False,
                  integration='verlet',myoutput=None, save_path=None):
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
    import gemmi
    from locscale.include.emmer.ndimage.map_tools import compute_real_space_correlation
    from locscale.include.emmer.pdb.pdb_to_map import pdb2map
    from locscale.include.emmer.ndimage.profile_tools import compute_radial_profile
    from locscale.include.emmer.pdb.pdb_utils import set_atomic_bfactors
    from locscale.include.emmer.pdb.modify_pdb import set_pdb_cell_based_on_gradient
    from locscale.preprocessing.pseudomodel_classes import Vector, add_Vector
    from locscale.utils.plot_tools import tab_print
    from tqdm import tqdm

    tabbed_print = tab_print(2)
    tprint = tabbed_print.tprint
    peak_bond_length_list = []
    map_values = []
    pseudomodel = model_initial.copy()
    gradient_magnitude = np.sqrt(gx**2+gy**2+gz**2)
    
    solver_properties = 'Solver started with the following properties: \n'+'\n Number of atoms = '+str(len(pseudomodel.list))+'\n Map potential: \n'+'\n g = '+str(g)+'\n Max gradient magnitude  = '+str(gradient_magnitude.max())+'\n Map value range  = '+str((emmap.min(),emmap.max()))+'\n Cap magnitude at  = '+str(capmagnitude_map)+'\n LJ Potential: \n'+'\n Equilibrium distance = '+str(min_dist_in_angst)+'\n Voxelsize, in A = '+str(voxelsize)+'\n LJ Factor = '+str(lj_factor)+'\n Epsilon = '+str(epsilon)+'\n Cap magnitude at  = '+str(capmagnitude_lj)+'\n Friction: \n'+ '\n Friction Coefficient = '+str(friction)+'\n Solver properties: \n'+'\n Total Iterations = '+str(total_iterations)+'\n Time step = '+str(dt)
    ## Convert solver_properties to a dictionary
    solver_properties_dict = {}
    for line in solver_properties.split('\n'):
        if '=' in line:
            key,value = line.split('=')
            solver_properties_dict[key.strip()] = value.strip()

    
        ## print the solver properties in a nice format
    
    print('='*50,file=myoutput)
    print('Solver started with the following properties: ',file=myoutput)
    for key,value in solver_properties_dict.items():
        print(key+' = '+value, file=myoutput)
    print('='*50,file=myoutput)
   
    if verbose:    
        print('# | Inter-atomic distance | Average map value | Average gradient acc | Average LJ potential', file=myoutput)        
    profiles_iterations = []
    cross_correlation = []
    for iter in tqdm(range(total_iterations),desc="Building Pseudo-atomic model"):
        
        neighborhood = get_neighborhood(pseudomodel.list,min_dist_in_angst/voxelsize)
        all_bond_lengths = np.array([d[0]*voxelsize for d in neighborhood.values()])
        bond_length_histogram = np.histogram(all_bond_lengths,bins=200)
        peak_bond_length = bond_length_histogram[1][bond_length_histogram[0].argmax()]
        peak_bond_length_list.append(peak_bond_length)
        gradient_list= []
        lj_list= []
        
        point_id = 0
        for atom in pseudomodel.list:            
            lj_neighbors = [pseudomodel.list[k] for k in neighborhood[point_id][1]]
            
            gradient_acceleration,map_value = get_acceleration_from_gradient(gx,gy,gz,emmap, g, point=atom, capmagnitude_map=capmagnitude_map)
            if len(lj_neighbors)==0:
                lj_potential_acceleration,lj_potential = Vector(np.array([0,0,0])),0
            else:
                lj_potential_acceleration,lj_potential = get_acceleration_from_lj_potential(atom, lj_neighbors, epsilon=1, min_dist_in_pixel=min_dist_in_angst/voxelsize,lj_factor=lj_factor,capmagnitude_lj=capmagnitude_lj)
            
            gradient_acceleration,lj_potential_acceleration = gradient_acceleration.scale(scale_map),lj_potential_acceleration.scale(scale_lj)
            acceleration = add_Vector(gradient_acceleration,lj_potential_acceleration)
            # add friction 
            atom.acceleration = add_Vector(acceleration, atom.velocity.scale(-friction))
            atom.map_value = map_value
            point_id += 1
            
            try:
                gradient_list.append(gradient_acceleration.magnitude())
                lj_list.append(lj_potential_acceleration.magnitude())
            except:
                gradient_list.append(0)
                lj_list.append(0)
            #lj_list
        
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

        if compute_map:
            pseudomodel.voxelsize = voxelsize
            pseudomodel.update_pdb_positions(voxelsize)
            gemmi_model=pseudomodel.convert_to_gemmi_model()
            
            emmap_shape = emmap.shape
            unitcell = gemmi.UnitCell(emmap_shape[0]*voxelsize,emmap_shape[1]*voxelsize,emmap_shape[2]*voxelsize,90,90,90)
            
            gemmi_structure = gemmi.Structure()
            gemmi_structure.add_model(gemmi_model)
            gemmi_structure.cell = unitcell
            gemmi_structure_bzero = set_atomic_bfactors(input_gemmi_st=gemmi_structure, b_iso=20)
            if save_path is not None:
                import os
                save_file = os.path.join(save_path,"pseudomodel_iteration_{}.pdb".format(iter))
                gemmi_structure_bzero.write_pdb(save_file)
                set_pdb_cell_based_on_gradient(gemmi_structure_bzero, emmap, apix=voxelsize, outpdb_name=os.path.join(save_path,"pseudomodel_iteration_{}.pdb".format(iter)))
            map_iteration = pdb2map(input_pdb=gemmi_structure,apix=voxelsize,size=emmap.shape, align_output=True)
            size = map_iteration.size
            cc = compute_real_space_correlation(map_iteration, emmap)
            cross_correlation.append(cc.max())
            profiles_iterations.append(compute_radial_profile(map_iteration))
            
    
            if verbose: 
                print(str(iter)+": #peak_bond_length = "+str(peak_bond_length_list[iter])\
                +": #map_value = "+str(map_values[iter])+": cc_max = "+str(cc.max()),file=myoutput)
    
        else:
            if verbose:
                gradient_acc_arr = np.array(gradient_list)
                lj_acc_arr = np.array(lj_list)
                print(str(iter)+
                      "\t | \t "+str(round(peak_bond_length_list[iter],2))+
                      "\t | \t "+str(map_values[iter])+
                      "\t | \t "+str(round(gradient_acc_arr.mean(),2))+
                      "\t | \t "+str(round(lj_acc_arr.mean(),2)),file=myoutput)    
            
    pseudomodel.voxelsize = voxelsize
    pseudomodel.update_pdb_positions(voxelsize)
    if returnPointsOnly == True:
        return pseudomodel    
    else:
        if compute_map:
            return pseudomodel, peak_bond_length_list, map_values, cross_correlation, profiles_iterations
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
    import random
    
    N = len(kicklist)
    points_array = np.array(points_array,dtype=float)
    for i in range(N):
        points_array[kicklist[i]]+=[random.uniform(-kick,kick),random.uniform(-kick,kick),random.uniform(-kick,kick)]
    return points_array    

        
def main_solver_kick(model_initial, min_dist_in_angst, voxelsize, total_iterations=99,returnPointsOnly=True,verbose=False, myoutput=None):
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
    from locscale.preprocessing.pseudomodel_classes import Model, Atom
    points_array = np.array([x.position.get() for x in model_initial.list])
    number_of_contacts = []
    if verbose:
        print(' Solver started with the following properties: \n'+
              '\n Number of atoms = '+str(len(points_array))+
              '\n Equilibrium distance = '+str(min_dist_in_angst)+
              '\n Voxelsize, in A = '+str(voxelsize)+
              '\n Total Iterations = '+str(total_iterations),file=myoutput)
              
    for i in range(total_iterations):
        neighbors = get_neighborhood(points_array,min_dist_in_angst,fromArray=True)
        kicklist = [x for x in neighbors.keys() if neighbors[x][0] <= min_dist_in_angst/voxelsize]
        points_array = find_and_kick(points_array,kicklist,kick=1)
        number_of_contacts.append(len(kicklist))
        if verbose: 
            print("Iteration number =  "+str(i)+": # Atoms less than eq. dist = "+str(len(kicklist)),file=myoutput)
        
        if sum(number_of_contacts[-3:]) == 0:
            break
        
    pseudomodel = Model([Atom(x) for x in points_array])
    pseudomodel.voxelsize = voxelsize
    pseudomodel.update_pdb_positions(voxelsize)
    if returnPointsOnly:
        return pseudomodel
    else:
        return pseudomodel, number_of_contacts


    