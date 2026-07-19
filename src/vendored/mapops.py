"""Map and cube utilities, vendored verbatim from LocScale.

Extracted with AST from locscale/include/emmer/ndimage/{map_utils,map_tools,filter}.py and
locscale/emmernet/emmernet_functions.py so the bundle is self-contained. The bodies are
unmodified; only intra-locscale imports were rewritten, since every referenced name now
lives in this module.
"""
import os

import numpy as np

# --- verbatim from locscale/include/emmer/ndimage/map_utils.py ---
def load_map(map_path, return_apix = True, verbose=False):
    import mrcfile
    pass  # (same module)
    emmap = mrcfile.open(map_path).data
    apix = average_voxel_size(mrcfile.open(map_path).voxel_size)
    
    if verbose:
        print("Loaded map from path: ", map_path)
        print("Voxel size: ", apix)
        print("Map shape: ", emmap.shape)
        
    if return_apix:
        return emmap, apix
    else:
        return emmap

# --- verbatim from locscale/include/emmer/ndimage/map_utils.py ---
def parse_input(input_map, allow_any_dims=True):
    '''
    Function to detect type of input and return a emmap numpy array

    Parameters
    ----------
    input_map : str or numpy array or Mrc object
    string type input should be path/to/emmap.mrc    

    Returns
    -------
    emmap : numpy.ndarray

    '''
    import os
    import mrcfile
    if isinstance(input_map, np.ndarray):
        if not allow_any_dims:
            if len(input_map.shape) == 3:
                return input_map
            else:
                print("You have not input a 3-D numpy array, which cannot be a EM-map")
                return None
        else:
            return input_map
            
    elif isinstance(input_map, str):
        if os.path.exists(input_map):
            emmap = mrcfile.open(input_map).data
            return emmap
        else:
            print("You have not entered a proper path, or the requested file does not exist!")
            return None

# --- verbatim from locscale/include/emmer/ndimage/map_utils.py ---
def save_as_mrc(map_data,output_filename, apix=None,origin=None,verbose=False, header=None):
    '''
    Function to save a numpy array containing volume, as a MRC file with proper header

    Parameters
    ----------
    map_data : numpy.ndarray
        Volume data showing the intensities of the EM Map at different points

    apix : float or any iterable
        In case voxelsize in x,y,z are all equal you can also just pass one parameter. 
    output_filename : str
        Path to save the MRC file. Example: 'path/to/map.mrc'
    origin: float or any iterable, optional
        In case origin index in x,y,z are all equal you can also just pass one parameter. 

    Returns
    -------
    Saves MRC .

    '''
    import mrcfile

    with mrcfile.new(output_filename,overwrite=True) as mrc:
        mrc.set_data(np.float32(map_data))
        
        if header is not None:
            mrc.set_extended_header(header)
        
        else:
            if apix is not None:
                #apix_list = [apix['x'], apix['y'], apix['z']]
                ## apix can be either a float or a list. If it's a single number, then the function convert_to_tuple will use it three times
                apix_tuple = convert_to_tuple(apix, num_dims=3)
                rec_array_apix = np.rec.array(apix_tuple, dtype=[('x','<f4'),('y','<f4'),('z','<f4')])
                mrc.voxel_size = rec_array_apix
            else:
                print("Please pass a voxelsize value either as a float or an iterable")
                return 0
            
            if origin is not None:    
                origin_tuple = convert_to_tuple(origin,num_dims=3)
            else:
                origin_tuple = convert_to_tuple(input_variable=0,num_dims=3)
            rec_array_origin = np.rec.array(origin_tuple, dtype=[('x','<f4'),('y','<f4'),('z','<f4')])
            mrc.header.origin = origin_tuple
            
        if verbose:
            print("Saving as MRC file format with following properties: ")
            print("File name: ", output_filename)
            print("Voxel size", mrc.voxel_size)
            print("Origin", mrc.header.origin)
            print("Shape", mrc.data.shape)
            
        
    mrc.close()

# --- verbatim from locscale/include/emmer/ndimage/map_utils.py ---
def extract_window(im, center, size):
    '''
    Extract a square window at a given location. 
    The center position of the window should be provided.

    Parameters
    ----------
    im : numpy.ndarray
        3D numpy array
    center : tuple, or list, or numpy.array (size=3)
        Position of the center of the window
    size : int, even
        Total window size (edge to edge) as an even number
        (In future could be modified to include different sized window 
        in different directions)
        

    Returns
    -------
    window : numpy.ndarray
        3D numpy array of shape (size x size x size)

    '''
    z,y,x = center
    window = im[z-size//2:z+size//2, y-size//2:y+size//2, x-size//2:x+size//2]
    return window

# --- verbatim from locscale/include/emmer/ndimage/map_utils.py ---
def binarise_map(*args, **kwargs):
    return binarize_map(*args, **kwargs)

# --- verbatim from locscale/include/emmer/ndimage/map_utils.py ---
def average_voxel_size(voxel_size_record):
    apix_x = voxel_size_record.x
    apix_y = voxel_size_record.y
    apix_z = voxel_size_record.z
    
    average_apix = (apix_x+apix_y+apix_z)/3
    
    return average_apix

# --- verbatim from locscale/include/emmer/ndimage/map_utils.py ---
def compute_FDR_confidenceMap_easy(em_map, apix, window_size, fdr=1, lowPassFilter_resolution=None,remove_temp_files=True, folder = None, use_default_noise_box=False):
    from .confidenceMapUtil.confidenceMapMain import calculateConfidenceMap
    pass  # (same module)
    import os, shutil, time
    
    if folder is None:
        current_cwd = os.getcwd()
    else:
        current_cwd = folder
    
    if not use_default_noise_box:
        noise_box_coords = detect_noise_boxes(em_map)
        print("Noise box coordinates detected: ", noise_box_coords)
    else:
        noise_box_coords = None
    timestamp =  str(time.time())
    temp_dir = current_cwd + '/fdr_output_temp_'+timestamp
    os.mkdir(temp_dir)
    os.chdir(temp_dir)
    confidenceMap,locFiltMap,locScaleMap,binMap,maskedMap = calculateConfidenceMap(
        em_map=em_map,apix=apix,noiseBox=noise_box_coords,testProc=None,ecdf=None,
        lowPassFilter_resolution=lowPassFilter_resolution,method=None, 
        window_size=window_size,windowSizeLocScale=None, locResMap=None,
        meanMap=None,varMap=None,fdr=fdr,modelMap=None,stepSize=None,mpi=None)
    
    fdr_threshold = np.min(maskedMap[np.nonzero(maskedMap)])
    
    os.chdir(current_cwd)
    if remove_temp_files:
        print("Clearing temporary files")
        shutil.rmtree(temp_dir)
    return confidenceMap, fdr_threshold

# --- verbatim from locscale/include/emmer/ndimage/map_utils.py ---
def resample_map(emmap, emmap_size_new=None, apix=None, apix_new=None, order=1, assert_shape=None):
    '''
    Function to resample an emmap in real space using linear interpolation 

    Parameters
    ----------
    emmap : numpy.ndimage
        
    emmap_size_new : tuple 
        
    apix : float
        
    apix_new : float
        

    Returns
    -------
    resampled_emmap

    '''
    from scipy.ndimage import zoom
    if emmap_size_new is None:
        if apix is not None and apix_new is not None:
            resample_factor = apix/apix_new
        else:
            raise UserWarning("Provide either (1) current pixel size and new pixel size or (2) new emmap size")
    
    else:
        try:
            resample_factor = emmap_size_new[0] / emmap.shape[0]
        except:
            raise UserWarning("Please provide proper input: emmap_size_new must be a tuple")
    
    if assert_shape is not None:
        if isinstance(assert_shape, int):
            nx = assert_shape
        if isinstance(assert_shape, tuple):
            nx = assert_shape[0]
        if isinstance(assert_shape, list):
            nx = assert_shape[0]
        assertion_factor = nx / (emmap.shape[0] * resample_factor)
        resample_factor *= assertion_factor

    resampled_image = zoom(emmap, resample_factor, order=order, grid_mode=False)


    
    return resampled_image

# --- verbatim from locscale/include/emmer/ndimage/map_tools.py ---
def detect_noise_boxes(emmap, num_windows=100):
    pass  # (same module)
    # find random centers
    emmap_shape = emmap.shape
    window_shape = int(emmap_shape[0] * 0.1) if emmap_shape[0] > 210 else 21
    emmap_shape = emmap.shape
    random_centers = get_random_center_voxels(window_shape, num_windows, emmap_shape)

    max_intensities_within_each_center = []
    for center in random_centers:
        window = extract_window(emmap, center, window_shape)
        max_intensities_within_each_center.append(np.max(window))
    
    index_of_center_with_least_max_intensity = np.argmin(max_intensities_within_each_center)
    center_with_least_max_intensity = random_centers[index_of_center_with_least_max_intensity]

    return center_with_least_max_intensity

# --- verbatim from locscale/include/emmer/ndimage/filter.py ---
def get_cosine_mask(mask,length_cosine_mask_1d):
    from scipy import signal
    cosine_window_1d = signal.windows.cosine(length_cosine_mask_1d)
    cosine_window_3d = window3D(cosine_window_1d)
    cosine_mask = signal.convolve(mask,cosine_window_3d,mode='same')
    cosine_mask = cosine_mask/cosine_mask.max()
    return cosine_mask

# --- verbatim from locscale/emmernet/emmernet_functions.py ---
def standardize_map(im):
    """ standardizes 3D density data

    Args:
        im (np.ndarray): 3D density data

    Returns:
        im (np.ndarray): standardized 3D density data
    """
    
    im = (im - im.mean()) / (10 * im.std())
    
    return im 

# --- verbatim from locscale/emmernet/emmernet_functions.py ---
def load_smoothened_mask(mask_path, mask_threshold=0.5, cosine_filter=3, verbose=False):
    pass  # (same module)
    pass  # (same module)
    
    mask, apix = load_map(mask_path, verbose=verbose)
    mask_binarize = (mask >= mask_threshold).astype(np.int_)
    mask_smooth = get_cosine_mask(mask_binarize, cosine_filter)
    mask_binarize = (mask_smooth >= mask_threshold).astype(np.int_)
    
    print("Mask threshold: {}".format(mask_threshold))
    print("Cosine filter: {}".format(cosine_filter))
    return mask_binarize, apix

# --- verbatim from locscale/emmernet/emmernet_functions.py ---
def extract_all_cube_centers(im_input, step_size, cube_size):
    '''
    Utility function to extract all cube centers from a 3D density map in a rolling window fashion
    
    '''
    im_shape = im_input.shape[0]
    length, width, height = im_input.shape

    # extract centers of all cubes in the 3D map based on the step size
    cubecenters = []
    for i in range(0, length, step_size):
        for j in range(0, width, step_size):
            for k in range(0, height, step_size):
                # i,j,k are corner of the cube 
                # we need to find the center of the cube
                center_k = k + cube_size//2
                center_j = j + cube_size//2
                center_i = i + cube_size//2

                # check if the center is within the map
                if center_k < length and center_j < width and center_i < height:
                    center_within_map = True
                else:
                    center_within_map = False
                
                # check if bounding box is within the map
                if k + cube_size < length and j + cube_size < width and i + cube_size < height:
                    bounding_box_within_map = True
                else:
                    bounding_box_within_map = False
                
                if center_within_map and bounding_box_within_map:
                    cubecenters.append((center_i, center_j, center_k))
                
                if center_within_map and not bounding_box_within_map:
                    # Check which dimension is out of bounds
                    if k + cube_size >= length:
                        diff  = k + cube_size - length
                        center_k = center_k - diff
                    if j + cube_size >= width:
                        diff  = j + cube_size - width
                        center_j = center_j - diff
                    if i + cube_size >= height:
                        diff  = i + cube_size - height
                        center_i = center_i - diff
                    cubecenters.append((center_i, center_j, center_k))
    
    return cubecenters

# --- verbatim from locscale/emmernet/emmernet_functions.py ---
def filter_cubecenters_by_mask(cubecenters, mask, cube_size, signal_to_noise_cubes, only_signal_cubes=False):
    '''
    Utility function to filter cube centers by a mask

    '''
    pass  # (same module)
    import random

    print("Initial number of cubes: {}".format(len(cubecenters)))
    filtered_cubecenters = []
    signal_cubes_centers = []
    noise_cubes_centers = []
    for center in cubecenters:
        cube = extract_window(mask, center=center, size=cube_size)
        if cube.sum() > 5:
            signal_cubes_centers.append(center)
        else:
            noise_cubes_centers.append(center)

    num_signal_cubes = len(signal_cubes_centers)
    num_noise_cubes = len(noise_cubes_centers)

    if only_signal_cubes:
        return signal_cubes_centers
    
    required_noise_cubes = int(num_signal_cubes / signal_to_noise_cubes)
    if num_noise_cubes < required_noise_cubes:
        print("Not enough noise cubes. Using all noise cubes")
        sampled_noise_cubes = noise_cubes_centers
        
    else:
        print(f"Using {required_noise_cubes} noise cubes out of {num_noise_cubes} noise cubes randomly")
        sampled_noise_cubes = random.sample(noise_cubes_centers, required_noise_cubes)
    print(f"num_signal_cubes: {num_signal_cubes}")
    print(f"num_noise_cubes: {len(sampled_noise_cubes)}")
    
    filtered_cubecenters = signal_cubes_centers + sampled_noise_cubes

    return filtered_cubecenters, signal_cubes_centers, sampled_noise_cubes

# --- verbatim from locscale/emmernet/emmernet_functions.py ---
def extract_cubes_from_cubecenters(emmap, cubecenters, cube_size):
    '''
    Utility function to extract all cubes from a 3D density map in a rolling window fashion
    
    '''
    pass  # (same module)
    import os
    import json
    # extract all cubes from the volume
    
    cubes = {}
    for i,center in enumerate(cubecenters):
        cube = extract_window(emmap, center=center, size=cube_size)
        cube = np.expand_dims(cube, axis=0)
    
        cubes[i] = {'cube': cube, 'center': center}
    
    # Extract the cubes array for input to the neural network
    cubes_array = np.array([cubes[i]['cube'] for i in range(len(cubes))])

    return cubes, cubes_array

# --- verbatim from locscale/emmernet/emmernet_functions.py ---
def get_cubes(emmap, step_size, cube_size, mask):
    
    # get cube centers
    cubecenters = extract_all_cube_centers(emmap, step_size, cube_size)
    filtered_signal_cubecenters = filter_cubecenters_by_mask(cubecenters, mask, cube_size, signal_to_noise_cubes=1, only_signal_cubes=True)
    
    cubes_dictionary, cubes_array = extract_cubes_from_cubecenters(emmap, filtered_signal_cubecenters, cube_size)
    
    return cubes_dictionary, cubes_array, filtered_signal_cubecenters

# --- verbatim from locscale/emmernet/emmernet_functions.py ---
def replace_cubes_in_dictionary(cubes_array, cubes_dictionary):
    cubes_dictionary_new = {}
    cubes_min = []
    cubes_max = []
    for i in range(len(cubes_array)):
        new_cube = cubes_array[i]
        if i > len(cubes_dictionary):
            continue
        cubes_dictionary_new[i] = {'cube': new_cube, 'center': cubes_dictionary[i]['center']}
        cubes_min.append(new_cube.min())
        cubes_max.append(new_cube.max())
    
    return cubes_dictionary_new

# --- verbatim from locscale/emmernet/emmernet_functions.py ---
def assemble_cubes(cubes_dictionary, im_shape, average=True):
    '''
    Utility function to assemble cubes into a 3D density map
    
    '''
    pass  # (same module)
    if isinstance(im_shape, int):
        imshape = (im_shape, im_shape, im_shape)
    else:
        imshape = im_shape
    
    im = np.zeros(imshape)
    average_map = np.zeros(imshape)
    for cubes in cubes_dictionary.values():
        center_ijk = cubes['center']
        ci, cj, ck = center_ijk

        cube = cubes['cube']
        if len(cube.shape) != 3:
            cube = cube.squeeze()
        ni, nj, nk = cube.shape

        im[ci-ni//2:ci+ni//2, cj-nj//2:cj+nj//2, ck-nk//2:ck+nk//2] += cube
        average_map[ci-ni//2:ci+ni//2, cj-nj//2:cj+nj//2, ck-nk//2:ck+nk//2] += 1
    
    if average:
        nonzero_indices = np.where(average_map != 0)
        im[nonzero_indices] /= average_map[nonzero_indices]

    return im


# --- verbatim from locscale/include/emmer/ndimage/map_utils.py ---
def binarize_map(emmap, threshold, return_type="int", threshold_type="gteq"):
    '''
    Function to binarize a map
    '''
    if threshold_type == "gteq":
        binary_map = emmap >= threshold
    elif threshold_type == "gt":
        binary_map = emmap > threshold
    elif threshold_type == "lteq":
        binary_map = emmap <= threshold
    elif threshold_type == "lt":
        binary_map = emmap < threshold
    else:
        print("Please provide a valid threshold_type")
        valid_threshold_types = ["gteq (>=)", "gt (>)", "lteq (<=)", "lt (<)"]
        raise ValueError(f"Invalid threshold_type provided {threshold_type}. Valid threshold_types are {valid_threshold_types}")

    if return_type == "int":
        binary_map = binary_map.astype(np.int_)
    elif return_type == "float":
        binary_map = binary_map.astype(np.float_)
    elif return_type == "bool":
        binary_map = binary_map.astype(bool)
    else:
        print("Please provide a valid return_type")
        valid_return_types = ["int", "float", "bool"]
        raise ValueError(f"Invalid return_type provided {return_type}. Valid return_types are {valid_return_types}")
    return binary_map


# --- verbatim from locscale/include/emmer/ndimage/map_utils.py ---
def convert_to_tuple(input_variable, num_dims=3):
    '''
    Convert any variable, or iterable into a tuple. If a scalar is input then a tuple is generated with same variable
    based on number of dimensions mentioned in num_dims

    Parameters
    ----------
    input_variable : any
        scalar, or any iterable
    num_dims : int, optional
        Length of tuple. The default is 3.
        
    Returns
    -------
    output_tuple : tuple

    '''
    
    if hasattr(input_variable, '__iter__'):
        if len(input_variable) == num_dims:
            output_tuple = tuple(input_variable)
            return output_tuple
        else:
            print("Input variable dimension {} doesn't match expected output dimension {}".format(len(input_variable), num_dims))
    else:
        output_list = [input_variable for temporary_index in range(num_dims)]
        output_tuple = tuple(output_list)
        return output_tuple


# --- verbatim from locscale/include/emmer/ndimage/map_tools.py ---
def get_random_center_voxels(window_shape, num_windows, emmap_shape):
    import random 
    random.seed(42)
    pass  # (same module)
    spherical_mask = get_spherical_mask(emmap_shape, emmap_shape[0]//2-window_shape)
    all_voxels_within_mask = np.asarray(np.where(spherical_mask == 1)).T.tolist()
    random_center_voxels = random.sample(all_voxels_within_mask, num_windows)
    return random_center_voxels


# --- verbatim from locscale/include/emmer/ndimage/filter.py ---
def window3D(w):
    # Convert a 1D filtering kernel to 3D
    # eg, window3D(numpy.hanning(5))
    
    L=w.shape[0]
    m1=np.outer(np.ravel(w), np.ravel(w))
    win1=np.tile(m1,np.hstack([L,1,1]))
    m2=np.outer(np.ravel(w),np.ones([1,L]))
    win2=np.tile(m2,np.hstack([L,1,1]))
    win2=np.transpose(win2,np.hstack([1,2,0]))
    win=np.multiply(win1,win2)
    return win


# --- verbatim from locscale/include/emmer/ndimage/filter.py ---
def get_spherical_mask(mask_shape, radius_index):
    n = mask_shape[0]
    z,y,x = np.ogrid[-n//2:n//2,-n//2:n//2,-n//2:n//2]
    mask = (x**2+y**2+z**2 <= radius_index**2).astype(int)
    return mask
