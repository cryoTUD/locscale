import numpy as np
import os
import keras
import time
from locscale.include.emmer.ndimage.map_utils import resample_map
#import gemmi


def compute_radial_profile_proper(vol, frequency_map):

    vol_fft = np.fft.rfftn(vol, norm="ortho");
    dim = vol_fft.shape;
    ps = np.real(np.abs(vol_fft));
    frequencies = np.fft.rfftfreq(dim[0]);
    #bins = np.digitize(frequency_map, frequencies);
    #bins = bins - 1;
    x, y, z = np.indices(ps.shape)
    radii = np.sqrt(x**2 + y**2 + z**2)
    radii = radii.astype(int)
    radial_profile = np.bincount(radii.ravel(), ps.ravel()) / np.bincount(radii.ravel())
    radial_profile = radial_profile[0:int(ps.shape[0]/2)+1]

    return radial_profile, frequencies;

def compute_scale_factors(em_profile, ref_profile):
    """Function to calculate the scale factors given two profiles.
    Returns:
        scale_factor (numpy array (1D)): The scale factors for the EM map
        bfactor (float): The local bfactor of the reference map 
        qfit (float): The local qfit of the reference map for bfactor calculation
    """
    
    ##############################################################################################
    # Stage 1: Calculate the scale factor
    ##############################################################################################

    np.seterr(divide='ignore', invalid='ignore')
    scale_factor = np.divide(np.abs(ref_profile), np.abs(em_profile))
    scale_factor[ ~ np.isfinite( scale_factor )] = 0; #handle division by zero    

    return scale_factor

def set_radial_profile(vol, scale_factors, frequencies, frequency_map, shape):
    vol_fft = np.fft.rfftn(np.copy(vol), norm='ortho');
    scaling_map = np.interp(frequency_map, frequencies, scale_factors);
    scaled_map_fft = scaling_map * vol_fft;
    scaled_map = np.real(np.fft.irfftn(scaled_map_fft, shape, norm='ortho'));

    return scaled_map, scaled_map_fft;

def get_central_scaled_pixel_vals_after_scaling_2(scaling_dictionary,verbose=False,process_name='LocScale'):

    """
    This function performs calls the scaling function in a rolling window fashion.
    Once the scaled cubes are calculated the central voxels in each cube is extracted
    into a list. This list is then returned. This function is compatible with
    both MPI and non-MPI environments.
    """
    from tqdm import tqdm
    from locscale.include.emmer.ndimage.map_tools import compute_real_space_correlation
    from locscale.include.emmer.ndimage.profile_tools import frequency_array, estimate_bfactor_standard
    from locscale.utils.math_tools import true_percent_probability
    from locscale.include.confidenceMapUtil import FDRutil
    import pickle
    from locscale.utils.math_tools import round_up_proper

    ###############################################################################
    # Stage 1: Initialize and collect variables
    ###############################################################################
    emmap = scaling_dictionary['emmap']
    masked_xyz_locs = scaling_dictionary['masked_xyz_locs']
    wn = scaling_dictionary['wn']
    apix = scaling_dictionary['apix']
    fsc_resolution = scaling_dictionary['fsc_resolution']
    processing_files_folder = scaling_dictionary['processing_files_folder']

    sharpened_vals = []
    qfit_voxels = []
    bfactor_voxels = []

    profiles_audit = {}

    temp_folder = processing_files_folder
    central_pix = round_up_proper(wn / 2.0)
    total = (masked_xyz_locs - wn / 2).shape[0]
    cnt = 1.0

    ###############################################################################
    # Stage 1a: Create a progress bar
    ###############################################################################

    mpi=False
    if process_name != 'LocScale':
        mpi=True
        from mpi4py import MPI

        comm = MPI.COMM_WORLD
        rank=comm.Get_rank()
        size=comm.Get_size()

        pbar = {}
        if rank == 0:
            description = "LocScale"
            pbar = tqdm(total=len(masked_xyz_locs)*size,desc=description)
    else:
        progress_bar=tqdm(total=len(masked_xyz_locs), desc=process_name)

    ###############################################################################
    # Stage 2: Perform the scaling in a rolling window fashion
    ###############################################################################
    frequency_map_window = FDRutil.calculate_frequency_map(np.zeros((wn, wn, wn)));
    freq = frequency_array(profile_size=25//2, apix=1)

    trained_model_path = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'trained_model_new_refined_models')
    reconstructed_model = keras.models.load_model(trained_model_path)
    #reconstructed_model = keras.models.load_model("/mnt/d/Modules/edited_locscale/locscale/locscale/utils/trained_model")

    for k, j, i in masked_xyz_locs - wn / 2:
        try:
            # k,j,i are indices of the corner voxel in each cube. Ensure it is rounded up to integer.
            k,j,i,wn = round_up_proper(k), round_up_proper(j), round_up_proper(i), round_up_proper(wn)

            #######################################################################
            # Stage 2a: Extract the cube from the EM map and model maps
            #######################################################################
            emmap_wn = emmap[k: k+wn, j: j+wn, i: i+ wn]


            #######################################################################
            # Stage 2b: Compute the radial profile of the two cubes
            #######################################################################
            em_profile, frequencies_map = compute_radial_profile_proper(emmap_wn, frequency_map_window)

            #########################
            # TO BE COMPLETED
            #########################

            # print(np.shape(em_profile))
            # print(em_profile)

            mod_profile = np.exp(reconstructed_model.predict(np.log(em_profile.reshape(1, 13, 1))))[0]

            # Checks scaling operation for 1% of all voxels.
            check_scaling=true_percent_probability(1)

            #######################################################################
            # Stage 2c: Compute the scale factors given the two radial profiles
            #######################################################################

            scale_factors = compute_scale_factors(em_profile, mod_profile)

            #######################################################################
            # Stage 2d: Get the scaled cube by applying the scale factors
            #######################################################################
            map_b_sharpened, map_b_sharpened_fft = set_radial_profile(emmap_wn, scale_factors, frequencies_map, frequency_map_window, emmap_wn.shape)

            #######################################################################
            # Stage 2e: For each cube, get the central voxel value of the scaled cube
            # and the bfactor information along with quality of fit
            #######################################################################

            bfactor, amplitude, qfit = estimate_bfactor_standard(
                        freq, em_profile, wilson_cutoff=10, fsc_cutoff=fsc_resolution, \
                        return_amplitude=True, return_fit_quality=True, standard_notation=True)
            if check_scaling:
                temp_dictionary = {}
                temp_dictionary['em_profile'] = em_profile
                temp_dictionary['ref_profile'] = mod_profile
                temp_dictionary['bfactor'] = bfactor
                temp_dictionary['amplitude'] = amplitude
                temp_dictionary['qfit'] = qfit
                temp_dictionary['scale_factors'] = scale_factors
                profiles_audit[tuple([k,j,i])] = temp_dictionary

            sharpened_vals.append(map_b_sharpened[central_pix, central_pix, central_pix])
            bfactor_voxels.append(bfactor)
            qfit_voxels.append(qfit)

        except Exception as e:

            #######################################################################
            # ERROR: If any error occurs, print the error and stop the operation
            #######################################################################
            print("Rogue voxel detected!  \n")
            print("Location (kji): {},{},{} \n".format(k,j,i))
            print("Skipping this voxel for calculation \n")
            k,j,i,wn = round_up_proper(k), round_up_proper(j), round_up_proper(i), round_up_proper(wn)

            emmap_wn = emmap[k: k+wn, j: j+wn, i: i+ wn]


            em_profile, frequencies_map = compute_radial_profile_proper(emmap_wn, frequency_map_window)


            print(em_profile)


            print(e)
            print(e.args)

            if mpi:
                print("Error occured at process: {}".format(rank))

            raise

        #### Progress bar update
        if mpi:
            if rank == 0:
                pbar.update(size)
        else:
            progress_bar.update(n=1)


    ###############################################################################
    # Stage 3: Save the processing files (the profile_audit file)
    ###############################################################################
    if mpi:
        if rank==0:
            import os
            pickle_file_output = os.path.join(temp_folder,"profiles_audit.pickle")
            with open(pickle_file_output,"wb") as audit:
                pickle.dump(profiles_audit, audit)
    else:

        import os

        pickle_file_output = os.path.join(temp_folder,"profiles_audit.pickle")
        with open(pickle_file_output,"wb") as audit:
            pickle.dump(profiles_audit, audit)

    ###############################################################################
    # Stage 4: Convert to numpy array and return the values
    ###############################################################################
    sharpened_vals_array = np.array(sharpened_vals, dtype=np.float32)
    bfactor_vals_array = np.array(bfactor_voxels, dtype=np.float32)
    qfit_vals_array = np.array(qfit_voxels, dtype=np.float32)


    return sharpened_vals_array , bfactor_vals_array, qfit_vals_array

def get_central_scaled_pixel_vals_after_scaling(scaling_dictionary,verbose=False,process_name='LocScale'):

    """ 
    This function performs calls the scaling function in a rolling window fashion. 
    Once the scaled cubes are calculated the central voxels in each cube is extracted
    into a list. This list is then returned. This function is compatible with 
    both MPI and non-MPI environments.
    """
    import os
    from tqdm import tqdm
    from locscale.include.emmer.ndimage.map_tools import compute_real_space_correlation
    from locscale.include.emmer.ndimage.profile_tools import frequency_array, estimate_bfactor_standard
    from locscale.utils.math_tools import true_percent_probability
    from locscale.include.confidenceMapUtil import FDRutil
    import pickle
    from locscale.utils.math_tools import round_up_proper
    import locscale

    
    
    ###############################################################################
    # Stage 1: Initialize and collect variables
    ###############################################################################
    emmap = scaling_dictionary['emmap']
    masked_xyz_locs = scaling_dictionary['masked_xyz_locs']
    wn = scaling_dictionary['wn']
    apix = scaling_dictionary['apix']
    fsc_resolution = scaling_dictionary['fsc_resolution']
    processing_files_folder = scaling_dictionary['processing_files_folder']
    
    sharpened_vals = []
    qfit_voxels = []
    bfactor_voxels = []
    
    profiles_audit = {}

    temp_folder = processing_files_folder
    central_pix = round_up_proper(wn / 2.0)
    total = (masked_xyz_locs - wn / 2).shape[0]
    cnt = 1.0

    ###############################################################################
    # Stage 1a: Create a progress bar 
    ###############################################################################
    
    progress_bar_assembly = tqdm(total=total, desc='Assembling radial profiles')
        
    ###############################################################################
    # Stage 2: Perform the scaling in a rolling window fashion
    ###############################################################################
    frequency_map_window = FDRutil.calculate_frequency_map(np.zeros((wn, wn, wn)));
    freq = frequency_array(profile_size=13, apix=1)

    trained_model_path = os.path.join(os.path.dirname(locscale.__file__),'utils','trained_model')
    reconstructed_model = keras.models.load_model(trained_model_path)

    raw_profiles = []
    
    
    for k, j, i in masked_xyz_locs - wn / 2:
        try:
           # print("Before rounding: {}".format((k,j,i)))
            # k,j,i are indices of the corner voxel in each cube. Ensure it is rounded up to integer.
            k,j,i,wn = round_up_proper(k), round_up_proper(j), round_up_proper(i), round_up_proper(wn)
          #  print("After rounding: {}".format((k,j,i)))
            
            #######################################################################
            # Stage 2a: Extract the cube from the EM map and model maps
            #######################################################################
            emmap_wn = emmap[k: k+wn, j: j+wn, i: i+ wn]
            

            #######################################################################
            # Stage 2b: Compute the radial profile of the two cubes        
            #######################################################################
            em_profile, frequencies_map = compute_radial_profile_proper(emmap_wn, frequency_map_window)
            raw_profiles.append(em_profile)
            #########################
            # TO BE COMPLETED
            #########################

            # print(np.shape(em_profile))
            # print(em_profile)

        except Exception as e:

                #######################################################################
                # ERROR: If any error occurs, print the error and stop the operation
                #######################################################################
                print("Rogue voxel detected!  \n")
                print("Location (kji): {},{},{} \n".format(k,j,i))
                print("Skipping this voxel for calculation \n")
                k,j,i,wn = round_up_proper(k), round_up_proper(j), round_up_proper(i), round_up_proper(wn)

                emmap_wn = emmap[k: k+wn, j: j+wn, i: i+ wn]


                em_profile, frequencies_map = compute_radial_profile_proper(emmap_wn, frequency_map_window)


                print(em_profile)


                print(e)
                print(e.args)

   

                raise
        
        progress_bar_assembly.update(1)

    em_profiles = np.array(raw_profiles)
    mod_profiles = np.exp(reconstructed_model.predict(np.log(em_profiles.reshape(np.shape(em_profiles)[0], 13, 1))))

    number = 0

    progress_bar_scaling = tqdm(total=total, desc='Scaling radial profiles')

    for k, j, i in masked_xyz_locs - wn / 2:
        try:
            k,j,i,wn = round_up_proper(k), round_up_proper(j), round_up_proper(i), round_up_proper(wn)
            emmap_wn = emmap[k: k+wn, j: j+wn, i: i+ wn]

            # Checks scaling operation for 1% of all voxels. 
            check_scaling=true_percent_probability(1) 
            
            #######################################################################
            # Stage 2c: Compute the scale factors given the two radial profiles
            #######################################################################
            mod_profile = mod_profiles[number]
            em_profile = raw_profiles[number]
            number = number + 1
            scale_factors = compute_scale_factors(em_profile, mod_profile)
            
            #######################################################################
            # Stage 2d: Get the scaled cube by applying the scale factors
            #######################################################################
            map_b_sharpened, map_b_sharpened_fft = set_radial_profile(emmap_wn, scale_factors, frequencies_map, frequency_map_window, emmap_wn.shape)
        
            #######################################################################
            # Stage 2e: For each cube, get the central voxel value of the scaled cube
            # and the bfactor information along with quality of fit
            #######################################################################
            
            bfactor, amplitude, qfit = estimate_bfactor_standard(
                        freq, em_profile, wilson_cutoff=10, fsc_cutoff=fsc_resolution, \
                        return_amplitude=True, return_fit_quality=True, standard_notation=True)
            if check_scaling:    
                temp_dictionary = {}
                temp_dictionary['em_profile'] = em_profile
                temp_dictionary['ref_profile'] = mod_profile
                temp_dictionary['bfactor'] = bfactor
                temp_dictionary['amplitude'] = amplitude
                temp_dictionary['freq'] = freq
                temp_dictionary['qfit'] = qfit
                temp_dictionary['scale_factors'] = scale_factors
                profiles_audit[tuple([k,j,i])] = temp_dictionary
                
            sharpened_vals.append(map_b_sharpened[central_pix, central_pix, central_pix])
            bfactor_voxels.append(bfactor)
            qfit_voxels.append(qfit)
            
        except Exception as e:

            #######################################################################
            # ERROR: If any error occurs, print the error and stop the operation
            #######################################################################
            print("Rogue voxel detected!  \n")
            print("Location (kji): {},{},{} \n".format(k,j,i))
            print("Skipping this voxel for calculation \n")
            k,j,i,wn = round_up_proper(k), round_up_proper(j), round_up_proper(i), round_up_proper(wn)
            
            emmap_wn = emmap[k: k+wn, j: j+wn, i: i+ wn]
            
        
            em_profile, frequencies_map = compute_radial_profile_proper(emmap_wn, frequency_map_window)
            
            
            print(em_profile)
            
                
            print(e)
            print(e.args)
            
            
            raise

        #### Progress bar update
        progress_bar_scaling.update(1)


    ###############################################################################
    # Stage 3: Save the processing files (the profile_audit file)
    ###############################################################################
    import os
    import datetime
    timestamp = datetime.datetime.now().strftime("%Y%m%d-%H%M%S")
    pickle_output_filename_with_timestamp = os.path.join(temp_folder, "profile_audit_{}.pickle".format(timestamp))
    with open(pickle_output_filename_with_timestamp,"wb") as audit:
        pickle.dump(profiles_audit, audit)

    ###############################################################################
    # Stage 4: Convert to numpy array and return the values
    ###############################################################################                                
    sharpened_vals_array = np.array(sharpened_vals, dtype=np.float32)
    bfactor_vals_array = np.array(bfactor_voxels, dtype=np.float32)
    qfit_vals_array = np.array(qfit_voxels, dtype=np.float32)
    
    results = {
        'sharpened_vals': sharpened_vals_array,
        'bfactor_vals': bfactor_vals_array,
        'qfit_vals': qfit_vals_array,
    }
    return results

def run_window_function_including_scaling(parsed_inputs_dict):
    """
    This is a function which performs high level data processing for Locscale

    """
    from locscale.utils.general import get_xyz_locs_and_indices_after_edge_cropping_and_masking
    from locscale.utils.general import save_list_as_map, put_scaled_voxels_back_in_original_volume_including_padding
    from locscale.utils.general import merge_sequence_of_sequences, split_sequence_evenly
    from joblib import Parallel, delayed
    ###############################################################################
    # Stage 1: Collect inputs from the dictionary
    ###############################################################################
    scaling_dictionary = parsed_inputs_dict
    
    ###############################################################################
    # Stage 2: Extract masked locations and indices from the mask
    ###############################################################################
    
    masked_xyz_locs, masked_indices, map_shape = get_xyz_locs_and_indices_after_edge_cropping_and_masking(
        scaling_dictionary['mask'], scaling_dictionary['wn'])

    masked_xyz_locs_split = split_sequence_evenly(masked_xyz_locs, scaling_dictionary['number_processes'])

    scaling_dictionary["masked_indices"] = masked_indices
    scaling_dictionary["map_shape"] = map_shape

    scaling_dictionary_split = {}
    for i in range(scaling_dictionary['number_processes']):
        scaling_dictionary_split[i] = scaling_dictionary.copy()
        scaling_dictionary_split[i]["masked_xyz_locs"] = masked_xyz_locs_split[i]
        scaling_dictionary_split[i]["use_joblib"] = True
    ###############################################################################
    # Stage 3: Run the window function to get sharpened values and bfactor information
    ###############################################################################
    # use joblib to run the window function in parallel
    if scaling_dictionary['number_processes'] > 1:
        results = Parallel(n_jobs=scaling_dictionary['number_processes'])(
            delayed(get_central_scaled_pixel_vals_after_scaling)(scaling_dictionary_split[i]) for i in range(scaling_dictionary['number_processes'])) 
    else:
        scaling_dictionary_split[0]["use_joblib"] = False
        results = [get_central_scaled_pixel_vals_after_scaling(scaling_dictionary_split[0])]
        
    #######################################################################
    # Stage 4: Merge the results from the parallel processing
    #######################################################################
    if scaling_dictionary['number_processes'] > 1:
        sharpened_vals = merge_sequence_of_sequences([results[i]['sharpened_vals'] for i in range(scaling_dictionary['number_processes'])])
        bfactor_vals = merge_sequence_of_sequences([results[i]['bfactor_vals'] for i in range(scaling_dictionary['number_processes'])])
        qfit_vals = merge_sequence_of_sequences([results[i]['qfit_vals'] for i in range(scaling_dictionary['number_processes'])])
    else:
        sharpened_vals = results[0]['sharpened_vals']
        bfactor_vals = results[0]['bfactor_vals']
        qfit_vals = results[0]['qfit_vals']
    
    ###############################################################################
    # Stage 4: Put the sharpened values back into the original volume
    ###############################################################################

    map_scaled = put_scaled_voxels_back_in_original_volume_including_padding(sharpened_vals, masked_indices, map_shape)
    
    ###############################################################################
    # Stage 5: Save processing files such as bfactor map and qfit maps 
    ###############################################################################
    
    bfactor_path = os.path.join(scaling_dictionary['processing_files_folder'], "bfactor_map.mrc")
    qfit_path = os.path.join(scaling_dictionary['processing_files_folder'], "qfit_map.mrc")
    save_list_as_map(bfactor_vals, masked_indices, map_shape, bfactor_path, scaling_dictionary['apix'])
    save_list_as_map(qfit_vals, masked_indices, map_shape, qfit_path, scaling_dictionary['apix'])

    ###############################################################################
    # Stage 6: Return the scaled map
    ###############################################################################

    # ADD post processing
    map_scaled = postprocess_map(map_scaled, scaling_dictionary['apix'])

    return map_scaled


def postprocess_map(emmap, apix):
    from locscale.include.emmer.ndimage.map_utils import resample_image
    from locscale.include.emmer.ndimage.map_utils import resample_map
    # resampling etc. and keep normilization
    emmap = resample_map(emmap, apix=1, apix_new=apix)
    #emmap = resample_image(emmap, apix=1, apix_new=apix)
    ## TBC
    return emmap

# locscale run_locscale -em /mnt/d/mapdata/0038_6gml/EMDBmaps/unsharpened_maps/emd_0038_unsharpened.map -res 3.2 -v -o /mnt/d/temp.mrc
