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

    #########################################################################
    # Import necessary modules
    #########################################################################
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

    #########################################################################
    # Stage 1: Check dependencies
    #########################################################################
    tabbed_print = tab_print(2)
    ## Check dependencies
    dependency_check = check_dependencies()
    if isinstance(dependency_check, list):
        print("The following dependencies are missing. The program may not work as expected. \n")
        print("\t".join(dependency_check))
    else:
        print("All dependencies are satisfied. \n")
    
    #########################################################################
    # Stage 2: Parse the inputs 
    # a) Prepare the emmap
    # b) Check axis orders of the maps 
    # c) Check if pseudo-model is required
    #########################################################################

    
    scale_using_theoretical_profile = not(args.ignore_profiles)  
    ##########################################################################
    ## scale_using_theoretical_profile is the flag used to determine 
    ## if scale factors are computed using the theoretical profile 
    ## which is required for the pseudo-atomic model routine. Although 
    ## pseudo-model routine is run automaticall if the pdb_path 
    ## is not provided, this flag can be used to override the theoretical 
    ## profile computation.
    ##########################################################################

    
    emmap_path, shift_vector = get_emmap_path_from_args(args)      
    xyz_emmap_path = check_axis_order(emmap_path)  
    xyz_emmap = mrcfile.open(xyz_emmap_path).data
    
    verbose = bool(args.verbose)
    fsc_resolution = float(args.ref_resolution)
    
    if args.apix is None:
        apix = average_voxel_size(mrcfile.open(emmap_path).voxel_size)  ## Assuming voxelsize is the same in all directions
    else:
        apix = float(args.apix)
    
    
    ###########################################################################
    # Stage 3: Prepare the mask
    ###########################################################################
    
    if verbose:
        print("."*80)
        print("Preparing mask \n")
    
    if args.mask is None:
        if args.verbose:
            tabbed_print.tprint("A mask path has not been provided. False Discovery Rate control (FDR) based confidence map will be calculated at 1% FDR \n")
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
    
    
    #############################################################################
    # Stage 4: Prepare the model-map
    # Here we check if the user has provided model map (.mrc format) or not. If the 
    # user has provided the model map, then we use it directly for computation. 
    # Else, we need to have a reference model and then simulate a model map from it. 
    # The reference model generation and simulation will be done in the 
    # preprocessing/pipeline module.
    #############################################################################

    ### Collect the required information for the reference model pipeline    
    if args.molecular_weight is not None:
        molecular_weight = float(args.molecular_weight)    
    else:
        molecular_weight = None

    if verbose:
        print("."*80)
        print("Preparing model map \n")

    ## If the user has not provided the model map and has not used the 
    ## no_reference option, then we need to run the get_modmap pipeline
    if args.model_map is None and not args.no_reference:  

        # Collect model map arguments and pass it to get_modmap pipeline
        
        pdb_path = args.model_coordinates
        
        ## Check if the user has provided the atomic model and set the
        ## scale_using_theoretical_profile to False if yes
        if pdb_path is not None:
            scale_using_theoretical_profile = False 
            ## If a PDB_path is provided, assume that it is an atomic model hence 
            ## we do not need to use the theoretical profile for scaling

            shift_coordinates(in_model_path=pdb_path, trans_matrix=shift_vector,
                                         out_model_path=pdb_path[:-4]+"_shifted.pdb")
            pdb_path = pdb_path[:-4]+"_shifted.pdb"

        ## Collect the rest of the arguments from argument parser
        add_blur = float(args.add_blur)
        skip_refine = args.skip_refine
        model_resolution = args.model_resolution
        ##### Following are arguments for pseudo-atomic model generation
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
        
        #############################################################################
        # Stage 4a: Pack all the collected arguments into a dictionary and pass it #
        #############################################################################
        
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
        
        #############################################################################
        # Stage 4b: Run the get_modmap pipeline                                 #
        #############################################################################
        modmap_path = get_modmap(modmap_args)
        xyz_modmap_path = check_axis_order(modmap_path, return_same_path=True)
        xyz_modmap = mrcfile.open(xyz_modmap_path).data
    
    ## If the user has provided the model map, then we use it directly for computation
    ## but not if the user has used the no_reference option

    elif args.model_map is not None and not args.no_reference:
        scale_using_theoretical_profile = False 
        ## If a model map is provide 
        ## we do not need to use the theoretical profile for scaling

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
        ## Running locscale without using any reference means that the bfactors of the 
        ## local window will be set to zero. This is not a recommended option.
        ## This option is only present due to testing purposes and may 
        ## be removed in the future.
        
        xyz_modmap = np.ones(xyz_emmap.shape)  ## only for code compliance
        set_local_bfactor = args.set_local_bfactor
    
    #############################################################################
    # Stage 5: Prepare other parameters for the locscale pipeline               
    #############################################################################
    if verbose:
        print("."*80)
        print("Preparing locscale parameters\n")

    ##############################################################################
    # Stage 5a: If window size is not given, then use a default size of 25A and 
    # calculate the window size in pixels based on this value
    ##############################################################################
    if args.window_size is None:   ## Use default window size of 25 A
        wn = round_up_to_even(25 / apix)
        if verbose:
            tabbed_print.tprint("Using a default window size of {} pixels, corresponding to approximately 25A".format(wn))
        
    elif args.window_size is not None:
        wn = round_up_to_even(int(args.window_size))
        if verbose:
            tabbed_print.tprint("Provided window size in pixels is {} corresponding to approximately {:.2f} Angstorm".format(wn, wn*apix))

    ##############################################################################
    # Stage 5b: If the locscale window extends out of the box then we need to
    # pad the input box to make it fit the window size
    ##############################################################################
    window_bleed_and_pad = check_for_window_bleeding(xyz_mask, wn)
    
    ## Collect the padded inputs if required
    if window_bleed_and_pad:
        pad_int_emmap = compute_padding_average(xyz_emmap, xyz_mask)
        pad_int_modmap = compute_padding_average(xyz_modmap, xyz_mask)
        map_shape = [(xyz_emmap.shape[0] + wn), (xyz_emmap.shape[1] + wn), (xyz_emmap.shape[2] + wn)]
        xyz_emmap = pad_or_crop_volume(xyz_emmap, map_shape, pad_int_emmap)
        xyz_modmap = pad_or_crop_volume(xyz_modmap, map_shape, pad_int_modmap)
        xyz_mask = pad_or_crop_volume(xyz_mask, map_shape, 0)


    ##############################################################################
    # Stage 5c: Get the parameters required to compute the scale factors
    # which maybe required for pseudo-atomic model routine when the 
    # scale_using_theoretical_profile is set to True
    ############################################################################## 

    ###############################################################################
    # Stage 6: Collect scale factor arguments and pack it into a dictionary
    ###############################################################################

    ###############################################################################
    # Definitions: 
    # wilson_cutoff: The resolution where Wilson regime starts (in Angstrom) 
    # high frequency cutoff: the resolution at which debye effects are negligible
    # FSC cutoff:  FSC resolution obtained using the two halfmaps thresholded at 0.143
    # boost_secondary_structure: Factor to boost radial profile of secondary structure 
    #                            in the debye regions
    # Nyquist cutoff: the resolution correponding to sampling frequency 
    # PS: all cutoff values are in Angstrom
    ###############################################################################
       
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
    
    ###############################################################################
    ## Add a dev_mode flag to override normal routine
    ###############################################################################
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
    
    ###############################################################################
    # Stage 6a: Pack into a dictionary
    ###############################################################################

    scale_factor_arguments = {}
    scale_factor_arguments['wilson'] = wilson_cutoff
    scale_factor_arguments['high_freq'] = high_frequency_cutoff
    scale_factor_arguments['fsc_cutoff'] = fsc_resolution
    scale_factor_arguments['nyquist'] = nyquist
    scale_factor_arguments['smooth'] = smooth_factor
    scale_factor_arguments['boost_secondary_structure'] = boost_secondary_structure
    scale_factor_arguments['no_reference'] = args.no_reference
    scale_factor_arguments['processing_files_folder'] = processing_files_folder
    if args.no_reference:
        scale_factor_arguments['set_local_bfactor'] = set_local_bfactor
    
    if verbose:
        print("Preparation completed. Now running LocScale!")
        print("."*80)
    
    
    #################################################################################
    # Stage 7: Pack everything into a dictionary and pass it to main function
    #################################################################################
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
    
    #################################################################################
    # Stage 8: Make some common sense checks and return 
    #################################################################################
    
    ## all maps should have same shape
    assert xyz_emmap.shape == xyz_modmap.shape == xyz_mask.shape, "The input maps and mask do not have the same shape"
    ## emmap and modmap should not be zeros
    assert abs(xyz_emmap.sum()) > 0 and abs(xyz_modmap.sum()) > 0, "Emmap and Modmap should not be zeros!"
    ## No element of the mask should be negative
    assert (xyz_mask>=0).any(), "Negative numbers found in mask"
    
    return parsed_inputs_dict
