import numpy as np

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
        

def load_smoothened_mask(mask_path, mask_threshold=0.5, cosine_filter=3, verbose=False):
    from locscale.include.emmer.ndimage.map_utils import load_map 
    from locscale.include.emmer.ndimage.filter import get_cosine_mask 
    
    mask, apix = load_map(mask_path, verbose=verbose)
    mask_binarize = (mask >= mask_threshold).astype(np.int_)
    mask_smooth = get_cosine_mask(mask_binarize, cosine_filter)
    mask_binarize = (mask_smooth >= mask_threshold).astype(np.int_)
    
    print("Mask threshold: {}".format(mask_threshold))
    print("Cosine filter: {}".format(cosine_filter))
    return mask_binarize, apix

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

def filter_cubecenters_by_mask(cubecenters, mask, cube_size, signal_to_noise_cubes, only_signal_cubes=False):
    '''
    Utility function to filter cube centers by a mask

    '''
    from locscale.include.emmer.ndimage.map_utils import extract_window
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


def extract_cubes_from_cubecenters(emmap, cubecenters, cube_size):
    '''
    Utility function to extract all cubes from a 3D density map in a rolling window fashion
    
    '''
    from locscale.include.emmer.ndimage.map_utils import extract_window , load_map, save_as_mrc
    import os
    import json
    # extract all cubes from the volume
    
    cubes = {}
    for i,center in enumerate(cubecenters):
        cube = extract_window(emmap, center=center, size=cube_size)
        cube = np.expand_dims(cube, axis=3)
    
        cubes[i] = {'cube': cube, 'center': center}
    
    # Extract the cubes array for input to the neural network
    cubes_array = np.array([cubes[i]['cube'] for i in range(len(cubes))])

    return cubes, cubes_array

def get_cubes(emmap, step_size, cube_size, mask):
    
    # get cube centers
    cubecenters = extract_all_cube_centers(emmap, step_size, cube_size)
    filtered_signal_cubecenters = filter_cubecenters_by_mask(cubecenters, mask, cube_size, signal_to_noise_cubes=1, only_signal_cubes=True)
    
    cubes_dictionary, cubes_array = extract_cubes_from_cubecenters(emmap, filtered_signal_cubecenters, cube_size)
    
    return cubes_dictionary, cubes_array, filtered_signal_cubecenters

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

    
def assemble_cubes(cubes_dictionary, im_shape, average=True):
    '''
    Utility function to assemble cubes into a 3D density map
    
    '''
    from locscale.include.emmer.ndimage.map_utils import extract_window
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

def show_signal_cubes(signal_cubes, im_shape, save_path, apix, input_shape):
    from locscale.include.emmer.ndimage.map_utils import save_as_mrc, resample_map
    emmap = np.zeros(im_shape)
    for cube in signal_cubes:
        emmap[cube[0], cube[1], cube[2]] = 1
    emmap_resampled = resample_map(emmap, apix=1, apix_new=apix, order=2, assert_shape=input_shape)
    save_as_mrc(emmap_resampled, save_path, apix=apix)

def calibrate_variance(variance_map):
    from locscale.emmernet.utils import load_calibrator
    
    calibrator = load_calibrator()
    input_map_shape = variance_map.shape
    calibrated_variance_map = calibrator.predict(variance_map.flatten()).reshape(input_map_shape)
    
    return calibrated_variance_map
            
        
def calculate_significance_map_from_emmernet_output(locscale_output_path, mean_prediction_path, var_prediction_path, n_samples=15):
    from locscale.include.emmer.ndimage.map_utils import load_map, save_as_mrc
    from locscale.emmernet.emmernet_functions import calibrate_variance
    from scipy.stats import norm
    import os
    import warnings
    # filter out warnings
    warnings.filterwarnings("ignore")
    # compute the z score map 
    locscale_map, apix = load_map(locscale_output_path)
    mean_prediction, apix = load_map(mean_prediction_path)
    var_prediction, apix = load_map(var_prediction_path)

    standard_deviation = np.sqrt(var_prediction)
    standard_error = standard_deviation / np.sqrt(n_samples)
    
    standard_error_calibrated = calibrate_variance(standard_error)

    z_score_map = (mean_prediction - locscale_map) / standard_error_calibrated

    # take care of nan values
    z_score_map[np.isnan(z_score_map)] = 100 # make it a large number so that is is very significant
    # convert the z score map to a p value map

    cdf_map = norm.cdf(z_score_map)
    # renormalize to -100 and 100
    pVDDT_map = cdf_map * 200 - 100
    
    
    output_folder = os.path.dirname(locscale_output_path)

    pVDDT_map_path = os.path.join(output_folder, "pVDDT.mrc")
    z_score_map_path = os.path.join(output_folder, "z_scores_calibrated.mrc")

    save_as_mrc(pVDDT_map, pVDDT_map_path, apix)
    save_as_mrc(z_score_map, z_score_map_path, apix)

    # Create a chimera script to visualize the map
    create_and_save_chimera_script(pVDDT_map_path, locscale_output_path, mean_prediction_path)

def create_and_save_chimera_script(pVDDT_map_path, locscale_output_path, mean_prediction_path):
    print("=========VISUALISATION INSTRUCTIONS=========")
    print("Open the feature enhanced map and pVDDT map in chimera")
    print("#1: Feature enhanced map and #2: pVDDT map")
    print("Copy the following commands:\n")
    threshold_command = "volume #1 level 0.12; volume #2 level 0.12;" 
    color_scheme_command = "color sample #1 map #2 palette -95,#0000ff:-80,#00ffff:0,#00ff00:80,#ffff00:95,#ff0000;" 
    print(threshold_command)
    print(color_scheme_command)
    print("============================================")
    