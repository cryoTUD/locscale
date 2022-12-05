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

       
    
    ##############################################################################
    # Stage 5b: If the locscale window extends out of the box then we need to
    # pad the input box to make it fit the window size
    ##############################################################################
    
    
    ##############################################################################
    # PreProcess the emmap 
    # Resample to 1A / pix
    # Normalise to mean 0 and std 0.1
    ##############################################################################
    
    preprocessed_emmap, preprocessed_mask = preprocess_map(xyz_emmap, xyz_mask, apix)
    wn = 26
    window_bleed_and_pad = check_for_window_bleeding(preprocessed_mask, wn)
    
    processing_files_folder = os.path.dirname(xyz_emmap_path)

    ## number of processes
    number_processes = args.number_processes
    
    if verbose:
        print("Preparation completed. Now running LocScale!")
        print("."*80)
    
    
    #################################################################################
    # Stage 7: Pack everything into a dictionary and pass it to main function
    #################################################################################
    parsed_inputs_dict = {}
    parsed_inputs_dict['emmap'] = preprocessed_emmap
    parsed_inputs_dict['mask'] = preprocessed_mask
    parsed_inputs_dict['wn'] = wn
    parsed_inputs_dict['apix'] = apix
 
    parsed_inputs_dict['win_bleed_pad'] = window_bleed_and_pad
    parsed_inputs_dict['fsc_resolution'] = fsc_resolution
    parsed_inputs_dict['emmap_path'] = xyz_emmap_path
    parsed_inputs_dict['mask_path'] = xyz_mask_path
    parsed_inputs_dict['processing_files_folder'] = processing_files_folder
    parsed_inputs_dict['number_processes'] = number_processes
    parsed_inputs_dict['verbose'] = verbose
    parsed_inputs_dict['original_map_shape'] = xyz_emmap.shape
    
    #################################################################################
    # Stage 8: Make some common sense checks and return 
    #################################################################################
    
    ## all maps should have same shape
    assert xyz_emmap.shape == xyz_mask.shape, "The input maps and mask do not have the same shape"
    ## emmap and modmap should not be zeros
    assert abs(xyz_emmap.sum()) > 0, "Emmap and Modmap should not be zeros!"
    ## No element of the mask should be negative
    assert (xyz_mask>=0).any(), "Negative numbers found in mask"
    
    return parsed_inputs_dict


    parsed_inputs_dict = {}
    
    return parsed_inputs_dict

def preprocess_map(emmap, mask, apix):
    from locscale.include.emmer.ndimage.map_utils import resample_image
    from locscale.include.emmer.ndimage.map_utils import resample_map
    from locscale.emmernet.run_emmernet import standardize_map

    emmap = resample_map(emmap, apix=apix, apix_new=1)
    mask = resample_map(mask, apix=apix, apix_new=1)

    #emmap = resample_image(emmap, apix=apix, apix_new=1)
    #mask = resample_image(mask, apix=apix, apix_new=1)

    emmap = standardize_map(emmap)
    
    return emmap, mask
