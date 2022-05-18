######################################################## IMPORTS ##################################################################

# external imports
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import mrcfile
import os
import math
from matplotlib.backends.backend_pdf import PdfPages

# internal imports
from locscale.include.emmer.ndimage.profile_tools import frequency_array
from locscale.include.emmer.ndimage.map_utils import average_voxel_size, resample_map

######################################################## FUNCTIONS ################################################################


def standardize_map(im):
    """ standardizes 3D density data

    Args:
        im (np.ndarray): 3D density data

    Returns:
        im (np.ndarray): standardized 3D density data
    """
    
    im = (im - im.mean()) / (10 * im.std())
    
    return im 


def minmax_normalize_map(im):
    """ normalizes 3D density data

    Args:
        im (np.ndarray): 3D density data

    Returns:
        im (np.ndarray): normalized 3D density data
    """
    
    im = (im - im.min()) / (im.max() - im.min())
    
    return im
        


def extract_single_cube_3D(cubes, cubecenters, cube_center, cube_size, im_input, i, verbose=True):
    """ extracts a single cube from the 3D density map

    Args:
        cubes (dict{int:np.ndarray}): contains all chunked cubes
        cubecenters (dict{int:dict{string:int}}): contains all chunked cube centers
        cube_center (dict{string:int}): the coordinate of the cube center in the coordinate system of the 3D map
        cube_size (int): cube shape in x, y, z direction
        im_input (np.ndarray): 3D density data
        i (int): specifies the cube and cube center number
        ax (axis): to plot the cube center positions on a visualized map slice
        verbose (bool): specifies whether or not to plot the cube center positions on the visualized map slice. Defaults to True.

    Returns:
        cubecenters (dict{int:dict{string:int}}): contains all chunked cube centers
        cube_center (dict{string:int}): the coordinate of the cube center in the coordinate system of the 3D map
        ax (axis): to plot the cube center positions on a visualized map slice
    """
    
    cube = im_input[int(cube_center['y'] - cube_size / 2) : int(cube_center['y'] + cube_size / 2),
                int(cube_center['x'] - cube_size / 2) : int(cube_center['x'] + cube_size / 2),
                int(cube_center['z'] - cube_size / 2) : int(cube_center['z'] + cube_size / 2),]   
    cubes[i+1] = cube
    cubecenters[i+1] = {'x': cube_center['x'], 'y': cube_center['y'], 'z': cube_center['z']}
    
    return cubes, cubecenters


def get_cubes(im_input, step_size, cube_size):
    """ extracts all cubes from a 3D density map

    Args:
        im_input (np.ndarray): 3D density data
        step_size (int): the step size in the x, y, z direction between chunked cubes
        cube_size (int): cube shape in x, y, z direction
        ax (axis): to plot the cube center positions on a visualized map slice

    Returns:
        cubes (np.ndarray): contains all chunked cubes
        cubecenters (np.ndarray): contains all chunked cube centers
        ax (axis): to plot the cube center positions on a visualized map slice
    """
    from tqdm import tqdm
    num_steps_x = int(np.floor(im_input.shape[0] / step_size))
    approximate_total_cubes = num_steps_x * num_steps_x * num_steps_x

    im_shape = im_input.shape[0]
    new_cube_center = {'x': cube_size / 2, 'y': cube_size / 2, "z": cube_size / 2} # first center
    cubes = {}
    cubecenters = {}
    last_round_z = last_round_y = x_end_reached = z_end_reached = y_end_reached = False
    i = 0
    
    progress_bar=tqdm(desc="Extracting cubes for prediction")
    while z_end_reached == False:
        while y_end_reached == False:
            while x_end_reached == False:
                # extract center
                progress_bar.update(1)
                cube_center = new_cube_center
                cubes, cubecenters = extract_single_cube_3D(cubes, cubecenters, cube_center, cube_size, im_input, i, True)
                
                # new center
                new_cube_center['x'] += step_size
                new_cube_end_x = new_cube_center['x'] + cube_size / 2
                i += 1
                
                # decide if / how to extract the next cube in the x-direction
                if new_cube_end_x <= im_shape:
                    pass
                elif im_shape < new_cube_end_x < im_shape + step_size:
                    new_cube_center['x'] -= (new_cube_end_x - im_shape)
                    cube_center = new_cube_center
                    cubes, cubecenters = extract_single_cube_3D(cubes, cubecenters, cube_center, cube_size, im_input, i, True)
                    i += 1
                    x_end_reached = True
                    if last_round_y:
                        y_end_reached = True
                    if last_round_z:
                        z_end_reached = True    
                elif new_cube_end_x >= im_shape + step_size:
                    x_end_reached = True
                    
            # reset x-position of the cube, set a next step in the y-direction
            new_cube_center['x'] = cube_size / 2
            new_cube_center['y'] += step_size
            x_end_reached = False
            
            # decide if / how to extract the next cube in the z-direction
            new_cube_end_y = new_cube_center['y'] + cube_size / 2
            if new_cube_end_y <= im_shape:
                pass
            elif im_shape < new_cube_end_y < im_shape + step_size:
                new_cube_center['y'] -= (new_cube_end_y - im_shape)
                last_round_y = True  
            elif new_cube_end_x >= im_shape + step_size:
                y_end_reached = True
            
        # reset x- and y-position of the cube, set a next step in the z-direction
        new_cube_center['x'] = cube_size / 2
        new_cube_center['y'] = cube_size / 2
        new_cube_center['z'] += step_size
        x_end_reached = False
        y_end_reached = False
        last_round_y = False
        
        # decide if / how to extract the next cube in the z-direction
        new_cube_end_z = new_cube_center['z'] + cube_size / 2
        if new_cube_end_z <= im_shape:
            pass
        elif im_shape < new_cube_end_z < im_shape + step_size:
            new_cube_center['z'] -= (new_cube_end_z - im_shape)
            last_round_z = True  
        elif new_cube_end_x >= im_shape + step_size:
            z_end_reached = True
    
    # dict to np array
    cubes = np.array(list(cubes.values()))
    cubecenters = np.array(list(cubecenters.values()))
    
    return cubes, cubecenters
    
   

def assemble_cubes(cubes, cubecenters, im_shape, cube_size):
    """ assembles chunked cubes back into a single map

    Args:
        cubes (np.ndarray): contains all chunked cubes
        cubecenters (np.ndarray): contains all chunked cube centers
        im_shape (int): map shape in x, y, z direction
        cube_size (int): cube shape in x, y, z direction

    Returns:
        im_reconstructed (np.ndarray): re-assembled 3D density map
        im_average_counter (np.ndarray): amount of density data averaging operations per map voxel
    """
    from tqdm import tqdm

    im_reconstructed, im_average_counter = np.zeros((im_shape, im_shape, im_shape)), np.zeros((im_shape, im_shape, im_shape))       
    
    for i in tqdm(np.arange(cubecenters.shape[0]), desc="Assembling cubes"):
        cube, cubecenter = cubes[i], cubecenters[i]
        top_left = {'y': int(cubecenter['y'] - cube_size / 2), 'x': int(cubecenter['x'] - cube_size / 2), 'z': int(cubecenter['z'] - cube_size / 2), }
        
     
        im_reconstructed[top_left['y']:top_left['y']+cube_size, top_left['x']:top_left['x']+cube_size, top_left['z']:top_left['z']+cube_size] += cube
        im_average_counter[top_left['y']:top_left['y']+cube_size, top_left['x']:top_left['x']+cube_size, top_left['z']:top_left['z']+cube_size] += 1

    
    im_reconstructed /= im_average_counter
    
    return im_reconstructed

