import numpy as np
import mrcfile
import locscale.include.emmer as emmer

def prepare_mask_and_maps_for_scaling(args):
    '''
    Parse the command line arguments and return inputs for computing local amplitude scaling 

    Parameters
    ----------
    args : Namespace

    Returns
    -------
    parsed_inputs_dict : dict
        Parsed inputs dictionary

    '''
    print("."*80)
    print("Preparing your inputs for LocScale")
    import os
    from locscale.preprocessing.pipeline import get_modmap
    from locscale.preprocessing.headers import run_FDR, check_axis_order
    from locscale.utils.math_tools import round_up_to_even, round_up_to_odd
    from locscale.utils.file_tools import get_emmap_path_from_args, check_dependencies
    from locscale.utils.general import get_spherical_mask, check_for_window_bleeding, compute_padding_average, pad_or_crop_volume
    from locscale.include.emmer.ndimage.map_tools import add_half_maps, compute_radial_profile_simple
    from locscale.include.emmer.ndimage.map_utils import average_voxel_size
    from locscale.include.emmer.ndimage.profile_tools import estimate_bfactor_through_pwlf, frequency_array, number_of_segments
    from locscale.include.emmer.pdb.pdb_tools import find_wilson_cutoff
    from locscale.include.emmer.pdb.pdb_utils import shift_coordinates
    from locscale.utils.plot_tools import tab_print
    
    tabbed_print = tab_print(2)
    ## Check dependencies
    dependency_check = check_dependencies()
    
    if isinstance(dependency_check, list):
        print("The following dependencies are missing. The program may not work as expected. \n")
        print("\t".join(dependency_check))
    else:
        print("All dependencies are satisfied. \n")
    
    scale_using_theoretical_profile = not(args.ignore_profiles)  ## Flag to determine whether to use theoretical profiles or not
    
    emmap_path, shift_vector = get_emmap_path_from_args(args)      
    xyz_emmap_path = check_axis_order(emmap_path)  
    
    xyz_emmap = mrcfile.open(xyz_emmap_path).data
    
    verbose = bool(args.verbose)
    
    fsc_resolution = float(args.ref_resolution)
    
    if args.apix is None:
        apix = average_voxel_size(mrcfile.open(emmap_path).voxel_size)  ## Assuming voxelsize is the same in all directions
    else:
        apix = float(args.apix)
    
    
    ## Get the mask path if provided. Calculate using FDR if a mask is not provided
    
    if verbose:
        print("."*80)
        print("Preparing mask \n")
    
    if args.mask is None:
        if args.verbose:
            tabbed_print.tprint("A mask path has not been provided. \
                 False Discovery Rate control (FDR) based confidence map will be calculated at 1% FDR \n")
        if args.fdr_window_size is None:   # if FDR window size is not set, take window size equal to 10% of emmap height
            fdr_window_size = round_up_to_even(xyz_emmap.shape[0] * 0.1)
            tabbed_print.tprint("FDR window size is not set. Using a default window size of {} \n".format(fdr_window_size))
        else:
            fdr_window_size = int(args.fdr_w)
        
        if args.fdr_filter is not None:
            filter_cutoff = float(args.fdr_filter)
            tabbed_print.tprint("A low pass filter value has been provided. \
                The EM-map will be low pass filtered to {:.2f} A \n".format(filter_cutoff))
        else:
            filter_cutoff = None
            
        mask_path = run_FDR(emmap_path=emmap_path, window_size = fdr_window_size, fdr=0.01, filter_cutoff=filter_cutoff)
        xyz_mask_path = check_axis_order(mask_path)
        
        if xyz_mask_path is not None:
            xyz_mask = (mrcfile.open(xyz_mask_path).data > 0.99).astype(np.int8)
        else:
            xyz_mask = get_spherical_mask(xyz_emmap.shape)
    else:
        mask_path = args.mask
        xyz_mask_path = check_axis_order(mask_path)
        xyz_mask = (mrcfile.open(xyz_mask_path).data > 0.99).astype(np.int8)
    
    
    ## Obtain the model map if provided if not determine from user input to generate pseudo-atomic model 

    
    if args.molecular_weight is not None:
        molecular_weight = float(args.molecular_weight)    
    else:
        molecular_weight = None

    if verbose:
        print("."*80)
        print("Preparing model map \n")

    if args.model_map is None and not args.no_reference:
        # Collect model map arguments and pass it to preprocessing pipeline
        #     
        pdb_path = args.model_coordinates
        if pdb_path is not None:
            scale_using_theoretical_profile = False ## If a PDB_path is provided, assume that it is an atomic model thus set this flag as False
            shift_coordinates(in_model_path=pdb_path, trans_matrix=shift_vector,
                                         out_model_path=pdb_path[:-4]+"_shifted.pdb")
            pdb_path = pdb_path[:-4]+"_shifted.pdb"
            
        add_blur = float(args.add_blur)
        skip_refine = args.skip_refine
        model_resolution = args.model_resolution
        ## Defaults for pseudo-atomic model 
        pseudomodel_method=args.pseudomodel_method
        pam_distance = float(args.distance)
        refmac_iter = int(args.refmac_iterations)
        refmac5_path = args.refmac5_path
        if pseudomodel_method == 'random' and args.total_iterations is None:
            pam_iteration = 100
        elif pseudomodel_method == 'gradient' and args.total_iterations is None:
            pam_iteration = 50
        elif args.total_iterations is not None:
            pam_iteration = int(args.total_iterations)
        build_ca_only = args.build_ca_only
        
        ## Get reference map using get_modmap_from_pseudomodel()
        ## Note that if a pdb_path is provided then the function 
        ## will use that instead of running pseudo-atomic model 
        ## routine. 
        
        modmap_args = {
            'emmap_path':xyz_emmap_path,
            'mask_path':xyz_mask_path,
            'pdb_path':pdb_path,
            'pseudomodel_method':pseudomodel_method,
            'pam_distance':pam_distance,
            'pam_iteration':pam_iteration,
            'fsc_resolution':fsc_resolution,
            'refmac_iter':refmac_iter,
            'add_blur':add_blur,
            'skip_refine':skip_refine,
            'model_resolution':model_resolution,
            'pg_symmetry':args.symmetry,
            'molecular_weight':molecular_weight,
            'build_ca_only':build_ca_only,
            'verbose':verbose,
            'refmac5_path':refmac5_path,
        }
        
        modmap_path = get_modmap(modmap_args)
        
        xyz_modmap_path = check_axis_order(modmap_path, return_same_path=True)
        xyz_modmap = mrcfile.open(xyz_modmap_path).data
    elif args.model_map is not None and not args.no_reference:
        scale_using_theoretical_profile = False ## If a model map is provided, assume that it is from an atomic model thus set this flag as False no matter what the user input 
        modmap_path = args.model_map
        model_resolution = args.model_resolution
        if model_resolution is not None:
            if verbose:
                tabbed_print.tprint("Performing low pass filter on the Model Map \
                    with a cutoff: {} based on user input".format(model_resolution))

            from locscale.include.emmer.ndimage.filter import low_pass_filter
            from locscale.include.emmer.ndimage.map_utils import save_as_mrc
            
            pseudo_map_unfiltered_data = mrcfile.open(modmap_path).data
            pseudo_map_filtered_data = low_pass_filter(im=pseudo_map_unfiltered_data, cutoff=model_resolution, apix=apix)
            
            filename = modmap_path[:-4]+"_filtered.mrc"
            save_as_mrc(map_data=pseudo_map_filtered_data, output_filename=filename, apix=apix)
            
            modmap_path = filename
        xyz_modmap_path = check_axis_order(modmap_path)
        xyz_modmap = mrcfile.open(xyz_modmap_path).data        
    else:
        print("Running locscale without using any reference")
        xyz_modmap = np.ones(xyz_emmap.shape)
    
    if verbose:
        print("."*80)
        print("Preparing locscale parameters\n")
    if args.window_size is None:   ## Use default window size of 25 A
        wn = round_up_to_even(25 / apix)
        if verbose:
            tabbed_print.tprint("Using a default window size of {} pixels, corresponding to approximately 25A".format(wn))
        
    elif args.window_size is not None:
        wn = round_up_to_even(int(args.window_size))
        if verbose:
            tabbed_print.tprint("Provided window size in pixels is {} corresponding to approximately {:.2f} Angstorm".format(wn, wn*apix))

    window_bleed_and_pad = check_for_window_bleeding(xyz_mask, wn)
    
    if window_bleed_and_pad:

        pad_int_emmap = compute_padding_average(xyz_emmap, xyz_mask)
        pad_int_modmap = compute_padding_average(xyz_modmap, xyz_mask)
        map_shape = [(xyz_emmap.shape[0] + wn), (xyz_emmap.shape[1] + wn), (xyz_emmap.shape[2] + wn)]
        xyz_emmap = pad_or_crop_volume(xyz_emmap, map_shape, pad_int_emmap)
        xyz_modmap = pad_or_crop_volume(xyz_modmap, map_shape, pad_int_modmap)
        xyz_mask = pad_or_crop_volume(xyz_mask, map_shape, 0)


    ## Next few lines of code characterizes radial profile of 
    ## input emmap : 
        ## wilson cutoff : threshold between guinier and wilson regimes in the radial profile
        ## high frequency cutoff : threshold above which to computing bfactor becomes valid  (for low resolution map, it's same as wilson cutoff)
        ## FSC cutoff : threshold above which amplitudes of signal becomes weaked compared to noise
        
    
    wilson_cutoff = find_wilson_cutoff(mask_path=xyz_mask_path, verbose=False)
    smooth_factor = args.smooth_factor
    boost_secondary_structure = args.boost_secondary_structure
    if fsc_resolution >= 6:
        high_frequency_cutoff = wilson_cutoff
        nyquist = (round(2*apix*10)+1)/10
        #fsc_cutoff = fsc_resolution
        bfactor_info = [0,np.array([0,0,0]),np.array([0,0,0])]
        pwlf_fit_quality = 0
    else:
        rp_emmap = compute_radial_profile_simple(xyz_emmap)
        freq = frequency_array(amplitudes=rp_emmap, apix=apix)
        num_segments = number_of_segments(fsc_resolution)
        bfactor, amp, (fit,z,slope) = estimate_bfactor_through_pwlf(freq=freq, amplitudes=rp_emmap, wilson_cutoff=wilson_cutoff, fsc_cutoff=fsc_resolution,num_segments=num_segments)
        nyquist = (round(2*apix*10)+1)/10
        #fsc_cutoff = fsc_resolution
        high_frequency_cutoff = 1/np.sqrt(z[-2])
        bfactor_info = [round(bfactor,2), 1/np.sqrt(z).round(2), np.array(slope).round(2)]  ## For information at end
        pwlf_fit_quality = fit.r_squared()
    
    dev_mode = args.dev_mode
    if dev_mode:
        print("DEV MODE: Scaling using theoretical profiles even if you have input an atomic model!")
        print("Before: scale_using_theoretical_profile=",scale_using_theoretical_profile)
        scale_using_theoretical_profile = None
        scale_using_theoretical_profile = True
        print("After: scale_using_theoretical_profile=",scale_using_theoretical_profile)        
    
    
        
    processing_files_folder = os.path.dirname(xyz_emmap_path)

    ## number of processes
    number_processes = args.number_processes
    
    scale_factor_arguments = {}
    scale_factor_arguments['wilson'] = wilson_cutoff
    scale_factor_arguments['high_freq'] = high_frequency_cutoff
    scale_factor_arguments['fsc_cutoff'] = fsc_resolution
    scale_factor_arguments['nyquist'] = nyquist
    scale_factor_arguments['smooth'] = smooth_factor
    scale_factor_arguments['boost_secondary_structure'] = boost_secondary_structure
    scale_factor_arguments['no_reference'] = args.no_reference
    scale_factor_arguments['processing_files_folder'] = processing_files_folder
    
    if verbose:
        print("Preparation completed. Now running LocScale!")
        print("."*80)
    
    
    parsed_inputs_dict = {}
    parsed_inputs_dict['emmap'] = xyz_emmap
    parsed_inputs_dict['modmap'] = xyz_modmap
    parsed_inputs_dict['mask'] = xyz_mask
    parsed_inputs_dict['wn'] = wn
    parsed_inputs_dict['apix'] = apix
 
    parsed_inputs_dict['use_theoretical'] = scale_using_theoretical_profile
    parsed_inputs_dict['scale_factor_args'] = scale_factor_arguments
    parsed_inputs_dict['verbose'] = verbose
    parsed_inputs_dict['win_bleed_pad'] = window_bleed_and_pad
    parsed_inputs_dict['bfactor_info'] = bfactor_info
    parsed_inputs_dict['fsc_resolution'] = fsc_resolution
    parsed_inputs_dict['PWLF_fit'] = pwlf_fit_quality
    parsed_inputs_dict['emmap_path'] = xyz_emmap_path
    parsed_inputs_dict['mask_path'] = xyz_mask_path
    parsed_inputs_dict['processing_files_folder'] = processing_files_folder
    parsed_inputs_dict['number_processes'] = number_processes
    
    
    ## all maps should have same shape
    assert xyz_emmap.shape == xyz_modmap.shape == xyz_mask.shape, "The input maps and mask do not have the same shape"
    ## emmap and modmap should not be zeros
    assert abs(xyz_emmap.sum()) > 0 and abs(xyz_modmap.sum()) > 0, "Emmap and Modmap should not be zeros!"
    ## No element of the mask should be negative
    assert (xyz_mask>=0).any(), "Negative numbers found in mask"
    
    return parsed_inputs_dict
