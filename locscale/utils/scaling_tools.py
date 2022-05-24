import numpy as np
import os
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

def compute_scale_factors(em_profile, ref_profile, apix, scale_factor_arguments, 
                          use_theoretical_profile=True, check_scaling=False):
    
    from locscale.include.emmer.ndimage.profile_tools import scale_profiles, merge_two_profiles, \
        add_deviations_to_reference_profile, frequency_array, estimate_bfactor_standard, get_theoretical_profile

    
    #print("checkScaling", check_scaling)
    #print("useTheoretical", use_theoretical_profile)
    if scale_factor_arguments['no_reference']:
        use_theoretical_profile = False
        from locscale.include.emmer.ndimage.profile_tools import estimate_bfactor_standard
        theoretical_profile_tuple = get_theoretical_profile(length=len(em_profile),apix=apix)
        freq = theoretical_profile_tuple[0]
        scaled_theoretical_amplitude = theoretical_profile_tuple[1]
        wilson_cutoff = scale_factor_arguments['wilson']
        fsc_cutoff = scale_factor_arguments['fsc_cutoff']
        b_sharpen, amplitude = estimate_bfactor_standard(freq, em_profile, wilson_cutoff, fsc_cutoff, return_amplitude=True)
        sharpening_profile = np.exp(0.25 * -1*b_sharpen * freq**2)
        scaled_reference_profile = em_profile * sharpening_profile
        fsc_filtered_reference_profile = merge_two_profiles(scaled_reference_profile, np.zeros(len(freq)), freq, smooth=10, d_cutoff=fsc_cutoff)
        
        deviated_reference_profile, exp_fit = add_deviations_to_reference_profile(freq, fsc_filtered_reference_profile, scaled_theoretical_amplitude, 
                                                                       wilson_cutoff=scale_factor_arguments['wilson'], 
                                                                       fsc_cutoff=scale_factor_arguments['nyquist'], 
                                                                       deviation_freq_start=scale_factor_arguments['wilson'], 
                                                                       deviation_freq_end=scale_factor_arguments['fsc_cutoff'], 
                                                                       magnify=scale_factor_arguments['boost_secondary_structure'])
        
        
        reference_profile_for_scaling = deviated_reference_profile
        if check_scaling:
            temporary_dictionary = {}
            temporary_dictionary['em_profile'] = em_profile
            temporary_dictionary['input_ref_profile'] = fsc_filtered_reference_profile
            temporary_dictionary['freq'] = freq
            temporary_dictionary['theoretical_amplitude'] = theoretical_profile_tuple[1]
            temporary_dictionary['scaled_theoretical_amplitude'] = scaled_theoretical_amplitude
            temporary_dictionary['scaled_reference_profile'] = scaled_reference_profile
            temporary_dictionary['fsc_filtered_reference_profile'] = fsc_filtered_reference_profile
            temporary_dictionary['deviated_reference_profile'] = deviated_reference_profile
            temporary_dictionary['exponential_fit'] = exp_fit
            temporary_dictionary['bfactor'] = b_sharpen
            temporary_dictionary['amplitude'] = amplitude
            temporary_dictionary['scaling_condition'] = scale_factor_arguments
    
    elif use_theoretical_profile:
        theoretical_profile_tuple = get_theoretical_profile(length=len(ref_profile),apix=apix)
        freq = theoretical_profile_tuple[0]
        
        num_atoms = ref_profile[0]
        mol_weight = num_atoms * 16  # daltons 
        wilson_cutoff_local = 1/(0.309 * np.power(mol_weight, -1/12))   ## From Amit Singer
        wilson_cutoff_local = np.clip(wilson_cutoff_local, scale_factor_arguments['fsc_cutoff']*1.5, scale_factor_arguments['wilson'])


        reference_profile_tuple = (freq, ref_profile)
        
        scaled_theoretical_tuple,(bfactor,amp, qfit) = scale_profiles(reference_profile_tuple, theoretical_profile_tuple,
                                                  wilson_cutoff=wilson_cutoff_local, fsc_cutoff=scale_factor_arguments['nyquist'], return_bfactor_properties=True)
        bfactor = -1 * bfactor  ## Standard notation
        
        scaled_theoretical_amplitude = scaled_theoretical_tuple[1]
        

        smooth = scale_factor_arguments['smooth']
        
        ## Using merge_profile
        scaled_reference_profile = merge_two_profiles(ref_profile,scaled_theoretical_amplitude,freq,smooth=smooth,d_cutoff=wilson_cutoff_local)
        
        ## Using deviations calculated directly from scaled theoretical profile
        deviations_begin = wilson_cutoff_local
        deviations_end = scale_factor_arguments['fsc_cutoff']
        magnify = scale_factor_arguments['boost_secondary_structure']
        
        deviated_reference_profile, exp_fit = add_deviations_to_reference_profile(freq, ref_profile, scaled_theoretical_amplitude, 
                                                                       wilson_cutoff=wilson_cutoff_local, 
                                                                       nyquist_cutoff=scale_factor_arguments['nyquist'], 

                                                                       deviation_freq_start=deviations_begin, 
                                                                       deviation_freq_end=deviations_end, 
                                                                       magnify=magnify)


        
        
        
        reference_profile_for_scaling = deviated_reference_profile
        if check_scaling:
            temporary_dictionary = {}
            temporary_dictionary['em_profile'] = em_profile
            temporary_dictionary['input_ref_profile'] = ref_profile
            temporary_dictionary['freq'] = freq
            temporary_dictionary['theoretical_amplitude'] = theoretical_profile_tuple[1]
            temporary_dictionary['scaled_theoretical_amplitude'] = scaled_theoretical_amplitude
            temporary_dictionary['scaled_reference_profile'] = scaled_reference_profile
            temporary_dictionary['deviated_reference_profile'] = deviated_reference_profile
            temporary_dictionary['exponential_fit'] = exp_fit
            temporary_dictionary['bfactor'] = bfactor
            temporary_dictionary['amplitude'] = amp
            temporary_dictionary['qfit'] = qfit
            temporary_dictionary['local_wilson'] = wilson_cutoff_local
            temporary_dictionary['deviations_begin'] = deviations_begin
            temporary_dictionary['deviations_end'] = deviations_end
            temporary_dictionary['magnify'] = magnify
            temporary_dictionary['scaling_condition'] = scale_factor_arguments
    else:
        freq = frequency_array(ref_profile, apix=apix)
        num_atoms = ref_profile[0]
        mol_weight = num_atoms * 16  # daltons 
        wilson_cutoff_local = 1/(0.309 * np.power(mol_weight, -1/12))   ## From Amit Singer
        wilson_cutoff_traditional = 10
        wilson_cutoff_local = np.clip(wilson_cutoff_local, scale_factor_arguments['fsc_cutoff']*1.5, scale_factor_arguments['wilson'])

        bfactor, amp, qfit = estimate_bfactor_standard(freq=freq, amplitude=ref_profile, wilson_cutoff=wilson_cutoff_traditional, 
                                                       fsc_cutoff=scale_factor_arguments['fsc_cutoff'], return_amplitude=True, return_fit_quality=True, standard_notation=True)
        
        reference_profile_for_scaling = ref_profile
        
        
    np.seterr(divide='ignore', invalid='ignore');
    
    scale_factor = np.divide(np.abs(reference_profile_for_scaling), np.abs(em_profile))
    scale_factor[ ~ np.isfinite( scale_factor )] = 0; #handle division by zero    
    
    if check_scaling and (use_theoretical_profile or scale_factor_arguments['no_reference']) :
        #print("checkScalingReport", check_scaling)
        temporary_dictionary['scale_factor'] = scale_factor
        return scale_factor, bfactor, qfit, temporary_dictionary
    else:
        return scale_factor, bfactor, qfit  

def set_radial_profile(vol, scale_factors, frequencies, frequency_map, shape):
    vol_fft = np.fft.rfftn(np.copy(vol), norm='ortho');
    scaling_map = np.interp(frequency_map, frequencies, scale_factors);
    scaled_map_fft = scaling_map * vol_fft;
    scaled_map = np.real(np.fft.irfftn(scaled_map_fft, shape, norm='ortho'));

    return scaled_map, scaled_map_fft;

def get_central_scaled_pixel_vals_after_scaling(emmap, modmap, masked_xyz_locs, wn, apix, use_theoretical_profile,scale_factor_arguments, verbose=False,f_cutoff=None, process_name='LocScale', audit=True):
    from tqdm import tqdm
    from locscale.include.emmer.ndimage.map_tools import compute_real_space_correlation
    from locscale.utils.math_tools import true_percent_probability
    from locscale.include.confidenceMapUtil import FDRutil
    import pickle
    from locscale.utils.math_tools import round_up_proper
    
    #if use_theoretical_profile:
    #    print("Using theoretical profiles for Local Scaling with the following parameters: \n")
    #    print(scale_factor_arguments)
    sharpened_vals = []
    qfit_voxels = []
    bfactor_voxels = []
    temp_folder = scale_factor_arguments["processing_files_folder"]
    
    central_pix = round_up_proper(wn / 2.0)
    total = (masked_xyz_locs - wn / 2).shape[0]
    cnt = 1.0
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
    
    frequency_map_window = FDRutil.calculate_frequency_map(np.zeros((wn, wn, wn)));

    if audit:
        profiles_audit = {}
    
    
    for k, j, i in masked_xyz_locs - wn / 2:
        
        try:
        
            k,j,i,wn = round_up_proper(k), round_up_proper(j), round_up_proper(i), round_up_proper(wn)
            
            emmap_wn = emmap[k: k+wn, j: j+wn, i: i+ wn]
            modmap_wn = modmap[k: k+wn, j: j+wn, i: i+ wn]
        
            em_profile, frequencies_map = compute_radial_profile_proper(emmap_wn, frequency_map_window);
            mod_profile, _ = compute_radial_profile_proper(modmap_wn, frequency_map_window);
                
            check_scaling=true_percent_probability(1) # Checks scaling operation for 1% of all voxels. 
                
            if check_scaling and use_theoretical_profile:
                scale_factors, bfactor, quality_fit, report = compute_scale_factors(em_profile, mod_profile,apix=apix,scale_factor_arguments=scale_factor_arguments, use_theoretical_profile=use_theoretical_profile,
        check_scaling=check_scaling)
                profiles_audit[(k,j,i)] = report
            else:
                scale_factors, bfactor, quality_fit = compute_scale_factors(em_profile, mod_profile,apix=apix, scale_factor_arguments=scale_factor_arguments, use_theoretical_profile=use_theoretical_profile,
        check_scaling=check_scaling)
                
                #map_b_sharpened = set_radial_profile(emmap_wn, scale_factors, radii)
            map_b_sharpened, map_b_sharpened_fft = set_radial_profile(emmap_wn, scale_factors, frequencies_map, frequency_map_window, emmap_wn.shape);
        
                #if verbose:
                #    if cnt%1000 == 0:
                #        print ('  {0} {1:.3} percent complete'.format(process_name, (cnt/total)*100))
                
        
            sharpened_vals.append(map_b_sharpened[central_pix, central_pix, central_pix])
            bfactor_voxels.append(bfactor)
            qfit_voxels.append(quality_fit)
            
        except Exception as e:
            print("Rogue voxel detected!  \n")
            print("Location (kji): {},{},{} \n".format(k,j,i))
            print("Skipping this voxel for calculation \n")
            k,j,i,wn = round_up_proper(k), round_up_proper(j), round_up_proper(i), round_up_proper(wn)
            
            emmap_wn = emmap[k: k+wn, j: j+wn, i: i+ wn]
            modmap_wn = modmap[k: k+wn, j: j+wn, i: i+ wn]
        
            em_profile, frequencies_map = compute_radial_profile_proper(emmap_wn, frequency_map_window);
            mod_profile, _ = compute_radial_profile_proper(modmap_wn, frequency_map_window);
            
            print(em_profile)
            print(mod_profile)
                
            print(e)
            print(e.args)
            
            if mpi:
                print("Error occured at process: {}".format(rank))
            
            raise
        
            
        if mpi:
            if rank == 0:
                pbar.update(size)
                
                    
        else:
            progress_bar.update(n=1)
            
    
        
    if mpi:
        if audit and use_theoretical_profile and rank==0:
            import os
            
            pickle_file_output = os.path.join(temp_folder,"profiles_audit.pickle")
            with open(pickle_file_output,"wb") as audit:
                pickle.dump(profiles_audit, audit)
    else:
        if audit and use_theoretical_profile:
            import os
            
            pickle_file_output = os.path.join(temp_folder,"profiles_audit.pickle")
            with open(pickle_file_output,"wb") as audit:
                pickle.dump(profiles_audit, audit)
                                
    sharpened_vals_array = np.array(sharpened_vals, dtype=np.float32)
    bfactor_vals_array = np.array(bfactor_voxels, dtype=np.float32)
    qfit_vals_array = np.array(qfit_voxels, dtype=np.float32)
    

    return sharpened_vals_array , bfactor_vals_array, qfit_vals_array

def run_window_function_including_scaling(parsed_inputs_dict):
    """
    >>> emmap, modmap, mask = setup_test_data()
    >>> scaled_vol = run_window_function_including_scaling(emmap,modmap,mask,wn=10,apix=1.0)
    >>> np.copy(EMNumPy.em2numpy(scaled_vol))[scaled_vol.get_xsize() / 2][scaled_vol.get_ysize() / 2]
    array([ 0.        ,  0.        ,  0.        ,  0.        ,  0.        ,
            0.        ,  0.        ,  0.        ,  0.        ,  0.        ,
            0.12524424,  0.15562208,  0.18547297,  0.24380369,  0.31203741,
            0.46546721,  0.47914436,  0.31334871,  0.28510684,  0.21345402,
            0.17892323,  0.        ,  0.        ,  0.        ,  0.        ,
            0.        ,  0.        ,  0.        ,  0.        ,  0.        ], dtype=float32)
    """
    mask = parsed_inputs_dict['mask']
    wn = parsed_inputs_dict['wn']
    emmap = parsed_inputs_dict['emmap']
    modmap = parsed_inputs_dict['modmap']
    use_theoretical_profile = parsed_inputs_dict['use_theoretical']
    scale_factor_arguments = parsed_inputs_dict['scale_factor_args']
    apix = parsed_inputs_dict['apix']
    verbose=parsed_inputs_dict['verbose']
    processing_files_folder=parsed_inputs_dict['processing_files_folder']
    from locscale.utils.general import get_xyz_locs_and_indices_after_edge_cropping_and_masking
    from locscale.utils.general import save_list_as_map, put_scaled_voxels_back_in_original_volume_including_padding
    
    masked_xyz_locs, masked_indices, map_shape = get_xyz_locs_and_indices_after_edge_cropping_and_masking(mask, wn)

    sharpened_vals, bfactor_vals, qfit_vals = get_central_scaled_pixel_vals_after_scaling(emmap, modmap, masked_xyz_locs, wn, apix, use_theoretical_profile,scale_factor_arguments=scale_factor_arguments,verbose=verbose)

    map_scaled = put_scaled_voxels_back_in_original_volume_including_padding(sharpened_vals, masked_indices, map_shape)
    
    

    ## Save temporary files 
    
    bfactor_path = os.path.join(processing_files_folder, "bfactor_map.mrc")
    qfit_path = os.path.join(processing_files_folder, "qfit_map.mrc")

        
    save_list_as_map(bfactor_vals, masked_indices, map_shape, bfactor_path, apix)
    save_list_as_map(qfit_vals, masked_indices, map_shape, qfit_path, apix)

    return map_scaled

def run_window_function_including_scaling_mpi(parsed_inputs_dict):
    """
    >>> emmap_name, modmap_name, mask_name = setup_test_data_to_files()
    >>> import subprocess
    >>> n = subprocess.call(mpi_cmd.split())
    >>> scaled_vol = get_image('scaled.mrc')
    >>> np.copy(EMNumPy.em2numpy(scaled_vol))[scaled_vol.get_xsize() / 2][scaled_vol.get_ysize() / 2]
    array([ 0.        ,  0.        ,  0.        ,  0.        ,  0.        ,
            0.        ,  0.        ,  0.        ,  0.        ,  0.        ,
            0.12524424,  0.15562208,  0.18547297,  0.24380369,  0.31203741,
            0.46546721,  0.47914436,  0.31334871,  0.28510684,  0.21345402,
            0.17892323,  0.        ,  0.        ,  0.        ,  0.        ,
            0.        ,  0.        ,  0.        ,  0.        ,  0.        ], dtype=float32)
    >>> n = [os.remove(each_file) for each_file in [emmap_name, modmap_name, mask_name, 'scaled.mrc']]
    """
    mask = parsed_inputs_dict['mask']
    wn = parsed_inputs_dict['wn']
    emmap = parsed_inputs_dict['emmap']
    modmap = parsed_inputs_dict['modmap']
    use_theoretical_profile = parsed_inputs_dict['use_theoretical']
    scale_factor_arguments = parsed_inputs_dict['scale_factor_args']
    apix = parsed_inputs_dict['apix']
    verbose=parsed_inputs_dict['verbose']
    processing_files_folder=parsed_inputs_dict['processing_files_folder']
    
    from mpi4py import MPI
    from locscale.utils.general import get_xyz_locs_and_indices_after_edge_cropping_and_masking
    from locscale.utils.general import save_list_as_map, merge_sequence_of_sequences, split_sequence_evenly, put_scaled_voxels_back_in_original_volume_including_padding
    
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    size = comm.Get_size()

    if rank == 0:
        masked_xyz_locs, masked_indices, map_shape = \
        get_xyz_locs_and_indices_after_edge_cropping_and_masking(mask, wn)

        zs, ys, xs = masked_xyz_locs.T
        zs = split_sequence_evenly(zs, size)
        ys = split_sequence_evenly(ys, size)
        xs = split_sequence_evenly(xs, size)
    else:
        zs = None
        ys = None
        xs = None

    zs = comm.scatter(zs, root=0)
    ys = comm.scatter(ys, root=0)
    xs = comm.scatter(xs, root=0)

    masked_xyz_locs = np.column_stack((zs, ys, xs))

    process_name = 'LocScale process {0} of {1}'.format(rank + 1, size)

    sharpened_vals, bfactor_vals, qfit_vals = get_central_scaled_pixel_vals_after_scaling(emmap, modmap, masked_xyz_locs, wn, apix,use_theoretical_profile=use_theoretical_profile, scale_factor_arguments=scale_factor_arguments,verbose=verbose,process_name=process_name)
    
    sharpened_vals = comm.gather(sharpened_vals, root=0)
    bfactor_vals = comm.gather(bfactor_vals, root=0)
    qfit_vals = comm.gather(qfit_vals, root=0)

    if rank == 0:
        sharpened_vals = merge_sequence_of_sequences(sharpened_vals)
        bfactor_vals = merge_sequence_of_sequences(bfactor_vals)
        qfit_vals = merge_sequence_of_sequences(qfit_vals)
        
        map_scaled = put_scaled_voxels_back_in_original_volume_including_padding(np.array(sharpened_vals),
        masked_indices, map_shape)
        
        print("Saving bfactor and qfist maps in here: {}".format(os.getcwd()))
            
        
        bfactor_path = os.path.join(processing_files_folder, "bfactor_map.mrc")
        qfit_path = os.path.join(processing_files_folder, "qfit_map.mrc")
            
        save_list_as_map(bfactor_vals, masked_indices, map_shape, bfactor_path, apix)
        save_list_as_map(qfit_vals, masked_indices, map_shape, qfit_path, apix)
        
    else:
        map_scaled = None

    comm.barrier()

    return map_scaled, rank




