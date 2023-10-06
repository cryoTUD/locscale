#
# Delft University of Technology (TU Delft) hereby disclaims all copyright interest in the program 'LocScale'
# written by the Author(s).
# Copyright (C) 2021 Alok Bharadwaj and Arjen J. Jakobi
# This software may be modified and distributed under the terms of the BSD license. 
# You should have received a copy of the BSD 3-clause license along with this program (see LICENSE file file for details).
# If not see https://opensource.org/license/bsd-3-clause/.
#

from json import load
import numpy as np
import locscale.include.emmer as emmer
from locscale.utils.plot_tools import tab_print
from locscale.preprocessing.headers import check_axis_order
from locscale.include.emmer.ndimage.map_utils import load_map
from locscale.utils.file_tools import RedirectStdoutToLogger, pretty_print_dictionary
tabbed_print = tab_print(2)
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
    args.logger.info("Preparing your inputs for LocScale")

    #########################################################################
    # If dev_mode is True, then the program will run in the development mode
    #########################################################################
    if args.dev_mode:
        return parse_inputs_from_processing_files(args)

    #########################################################################
    # Import necessary modules
    #########################################################################
    import os
    from locscale.utils.file_tools import get_emmap_path_from_args, check_dependencies   

    #########################################################################
    # Stage 1: Check dependencies
    #########################################################################
    
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
    parsed_inputs = vars(args)

    # Prepare the emmap
    args.logger.info("1) Parsing the unsharpened map inputs")
    parsed_inputs["unsharpened_emmap_path"], parsed_inputs["shift_vector"]  = get_emmap_path_from_args(args)
    # Check axis orders of the maps
    parsed_inputs["xyz_emmap_path"] = check_axis_order(parsed_inputs["unsharpened_emmap_path"])
    with RedirectStdoutToLogger(args.logger, wait_message="Loading inputs"):
        parsed_inputs["xyz_emmap"], apix_from_file = load_map(parsed_inputs["xyz_emmap_path"], verbose=True)
    
    parsed_inputs["apix"] = parsed_inputs["apix"] if parsed_inputs["apix"] else apix_from_file
    
    ##########################################################################
    ## use_theoretical_profile is the flag used to determine 
    ## if scale factors are computed using the theoretical profile 
    ## which is required for the pseudo-atomic model routine. Although 
    ## pseudo-model routine is run automaticall if the pdb_path 
    ## is not provided, this flag can be used to override the theoretical 
    ## profile computation.
    ##########################################################################
           
    ###########################################################################
    # Stage 3: Prepare the mask
    ###########################################################################
    
    if parsed_inputs["verbose"]:
        print("."*80)
        print("Preparing mask \n")
    args.logger.info("2) Preparing mask")
    parsed_inputs["xyz_mask"], parsed_inputs["xyz_mask_path"], parsed_inputs["mask_path_raw"] = prepare_mask_from_inputs(parsed_inputs)
    
    
    #############################################################################
    # Stage 4: Prepare the model-map
    # Here we check if the user has provided model map (.mrc format) or not. If the 
    # user has provided the model map, then we use it directly for computation. 
    # Else, we need to have a reference model and then simulate a model map from it. 
    # The reference model generation and simulation will be done in the 
    # preprocessing/pipeline module.
    #############################################################################

    if parsed_inputs["no_reference"]:
        print("Running locscale without using any reference")
        args.logger.info("3) Running locscale without using any reference")
        ## Running locscale without using any reference means that the bfactors of the 
        ## local window will be set to zero. This is not a recommended option.
        ## This option is only present due to testing purposes and may 
        ## be removed in the future.
        
        parsed_inputs["xyz_modmap"] = np.ones(parsed_inputs["xyz_emmap"].shape)  ## only for code compliance
        parsed_inputs["xyz_modmap_path"] = None

    else:
        print("."*80)
        print("Preparing model map \n")
        args.logger.info("3) Preparing model map")
        parsed_inputs["xyz_modmap"], parsed_inputs["xyz_modmap_path"] = get_modmap_from_inputs(parsed_inputs)


    #############################################################################
    # Stage 5: Prepare other parameters for the locscale pipeline               
    #############################################################################
    if parsed_inputs["verbose"]:
        print("."*80)
        print("Preparing locscale parameters\n")

    ##############################################################################
    # Stage 5a: If window size is not given, then use a default size of 25A and 
    # calculate the window size in pixels based on this value
    ##############################################################################
    
    parsed_inputs["wn"] = set_window_size(parsed_inputs["window_size"], parsed_inputs["apix"], parsed_inputs["verbose"])
    parsed_inputs["processing_files_folder"] = os.path.dirname(parsed_inputs["xyz_emmap_path"])

    ##############################################################################
    # Stage 5b: If the locscale window extends out of the box then we need to
    # pad the input box to make it fit the window size
    ##############################################################################
    parsed_inputs = pad_input_maps_if_required(parsed_inputs)
    
    ##############################################################################
    # Stage 5c: Get the parameters required to compute the scale factors
    # which maybe required for pseudo-atomic model routine when the 
    # use_theoretical_profile is set to True
    ############################################################################## 
    parsed_inputs["scale_factor_args"] = get_scale_factor_arguments(parsed_inputs)
    
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
       
    if parsed_inputs["verbose"]:
        print("Preparation completed. Now running LocScale!")
        parsed_inputs["logger"].info("Preparation completed. Now running LocScale!")
        print("."*80)
    
    
    #################################################################################
    # Stage 7: Pack everything into a dictionary and pass it to main function
    #################################################################################
    parsed_inputs_dict = {}
    parsed_inputs_dict['emmap'] = parsed_inputs["xyz_emmap"]
    parsed_inputs_dict['modmap'] = parsed_inputs["xyz_modmap"]
    parsed_inputs_dict['mask'] = parsed_inputs["xyz_mask"]
    parsed_inputs_dict['wn'] = parsed_inputs["wn"]
    parsed_inputs_dict['apix'] = parsed_inputs["apix"]
    parsed_inputs_dict['use_theoretical_profile'] = parsed_inputs["use_theoretical_profile"]
    parsed_inputs_dict['scale_factor_arguments'] = parsed_inputs["scale_factor_args"]
    parsed_inputs_dict['verbose'] = parsed_inputs["verbose"]
    parsed_inputs_dict['win_bleed_pad'] = parsed_inputs["window_bleed_and_pad"]
    parsed_inputs_dict['bfactor_info'] = parsed_inputs["scale_factor_args"]["bfactor_info"]
    parsed_inputs_dict['fsc_resolution'] = parsed_inputs["ref_resolution"]
    parsed_inputs_dict['PWLF_fit'] = parsed_inputs["scale_factor_args"]["pwlf_fit_quality"]
    parsed_inputs_dict['emmap_path'] = parsed_inputs["xyz_emmap_path"]
    parsed_inputs_dict['mask_path'] = parsed_inputs["xyz_mask_path"]
    parsed_inputs_dict['processing_files_folder'] = parsed_inputs["processing_files_folder"]
    parsed_inputs_dict['number_processes'] = parsed_inputs["number_processes"]
    parsed_inputs_dict['complete_model'] = parsed_inputs["complete_model"]
    parsed_inputs_dict['original_map_shape'] = parsed_inputs["original_map_shape"]
    
    
    parsed_inputs["logger"].info("Parsed inputs dictionary: \n{}".format(pretty_print_dictionary(parsed_inputs_dict)))
    
    #################################################################################
    # Stage 8: Make some common sense checks and return 
    #################################################################################
    
    ## all maps should have same shape
    assert parsed_inputs_dict["emmap"].shape == parsed_inputs_dict["modmap"].shape == parsed_inputs_dict['mask'].shape, "The input maps and mask do not have the same shape"
    ## emmap and modmap should not be zeros
    assert abs(parsed_inputs_dict["emmap"].sum()) > 0 and abs(parsed_inputs_dict["modmap"].sum()) > 0, "Emmap and Modmap should not be zeros!"
    ## No element of the mask should be negative
    assert (parsed_inputs_dict['mask']>=0).any(), "Negative numbers found in mask"
    
    # # Dump the parsed inputs to a pickle file in the input folder
    # import pickle
    # with open(os.path.join(parsed_inputs_dict['processing_files_folder'], 'parsed_inputs.pickle'), 'wb') as f:
    #     pickle.dump(parsed_inputs_dict, f)
        
    return parsed_inputs_dict

def parse_inputs_from_processing_files(args):
    from locscale.include.emmer.ndimage.map_utils import load_map, save_as_mrc
    from locscale.utils.file_tools import get_emmap_path_from_args
    from locscale.preprocessing.headers import check_axis_order
    import os
    import pickle

    original_emmap, _ = get_emmap_path_from_args(args)
    original_emmap = load_map(original_emmap)[0]
    original_map_shape = original_emmap.shape

    parsed_inputs = vars(args)
    processing_files_folder = parsed_inputs["output_processing_files"]
    parsed_inputs['processing_files_folder'] = processing_files_folder

    intermediate_output_pickle = os.path.join(processing_files_folder, "intermediate_outputs.pickle")
    assert os.path.exists(intermediate_output_pickle), "Intermediate output pickle file not found! Do not use dev_mode!"
    with open(intermediate_output_pickle, "rb") as f:
        intermediate_outputs = pickle.load(f)

    parsed_inputs["unsharpened_emmap_path"], parsed_inputs["shift_vector"]  = get_emmap_path_from_args(args)

    parsed_inputs["xyz_emmap_path"] = intermediate_outputs["emmap_path"]
    parsed_inputs["xyz_mask_path"] = intermediate_outputs["mask_path"]
    parsed_inputs["mask_path_raw"] = intermediate_outputs["mask_path"]
    parsed_inputs["pseudomodel_modmap"] = intermediate_outputs["pseudomodel_modmap"]
    parsed_inputs["difference_mask_path"] = intermediate_outputs["difference_mask_path"]

    xyz_emmap, apix = load_map(parsed_inputs["xyz_emmap_path"])
    xyz_mask, _ = load_map(parsed_inputs["xyz_mask_path"])
    xyz_modmap, _ = load_map(parsed_inputs["pseudomodel_modmap"])

    
    parsed_inputs["apix"] = apix
    parsed_inputs["original_map_shape"] = original_map_shape
    parsed_inputs["xyz_emmap"] = xyz_emmap
    parsed_inputs["xyz_mask"] = xyz_mask
    parsed_inputs["xyz_modmap"] = xyz_modmap
    
    parsed_inputs["pseudomodel_required"] = True if parsed_inputs["model_coordinates"] is None else False
    parsed_inputs["use_theoretical_profile"] = parsed_inputs["pseudomodel_required"]
    parsed_inputs["wn"] = set_window_size(parsed_inputs["window_size"], apix, parsed_inputs["verbose"])
    
    parsed_inputs = pad_input_maps_if_required(parsed_inputs)

    parsed_inputs["scale_factor_args"] = get_scale_factor_arguments(parsed_inputs)

    parsed_inputs_dict = {}
    parsed_inputs_dict['emmap'] = xyz_emmap
    parsed_inputs_dict['modmap'] = xyz_modmap
    parsed_inputs_dict['mask'] = xyz_mask
    parsed_inputs_dict['wn'] = parsed_inputs["wn"]
    parsed_inputs_dict['apix'] = apix
    parsed_inputs_dict['use_theoretical_profile'] = parsed_inputs["use_theoretical_profile"]
    parsed_inputs_dict['scale_factor_arguments'] = parsed_inputs["scale_factor_args"]
    parsed_inputs_dict['verbose'] = parsed_inputs["verbose"]
    parsed_inputs_dict['win_bleed_pad'] = parsed_inputs["window_bleed_and_pad"]
    parsed_inputs_dict['bfactor_info'] = parsed_inputs["scale_factor_args"]["bfactor_info"]
    parsed_inputs_dict['fsc_resolution'] = parsed_inputs["ref_resolution"]
    parsed_inputs_dict['PWLF_fit'] = parsed_inputs["scale_factor_args"]["pwlf_fit_quality"]
    parsed_inputs_dict['emmap_path'] = parsed_inputs["xyz_emmap_path"]
    parsed_inputs_dict['mask_path'] = parsed_inputs["xyz_mask_path"]
    parsed_inputs_dict['processing_files_folder'] = parsed_inputs["output_processing_files"]
    parsed_inputs_dict['number_processes'] = parsed_inputs["number_processes"]
    parsed_inputs_dict['complete_model'] = parsed_inputs["complete_model"]
    parsed_inputs_dict['original_map_shape'] = parsed_inputs["original_map_shape"]
    

    assert parsed_inputs_dict["emmap"].shape == parsed_inputs_dict["modmap"].shape == parsed_inputs_dict['mask'].shape, "The input maps and mask do not have the same shape"
    ## emmap and modmap should not be zeros
    assert abs(parsed_inputs_dict["emmap"].sum()) > 0 and abs(parsed_inputs_dict["modmap"].sum()) > 0, "Emmap and Modmap should not be zeros!"
    ## No element of the mask should be negative
    assert (parsed_inputs_dict['mask']>=0).any(), "Negative numbers found in mask"

    return parsed_inputs_dict

def choose_fdr_filter_based_on_inputs(parsed_inputs):
    if parsed_inputs["fdr_filter"] is not None:
        filter_cutoff = float(parsed_inputs["fdr_filter"])
    else:
        if parsed_inputs["run_type"] == "feature_enhance": # "feature_enhance" is a key passed when running emmernet prediction
            filter_cutoff = 5 # default filter cutoff for estimating mask for emmernet prediction
        else:
            if parsed_inputs["modality"] == "predict_model_map":
                filter_cutoff = 5
            else:
                filter_cutoff = None
                
    parsed_inputs["logger"].info("FDR filter cutoff is set to {} \n".format(str(filter_cutoff)))
    
    return filter_cutoff

def prepare_inputs_for_FDR(parsed_inputs):
    from locscale.utils.math_tools import round_up_to_even
    
    # Set fdr_window_size
    if parsed_inputs["fdr_window_size"] is None:   # if FDR window size is not set, take window size equal to 10% of emmap height
        fdr_window_size = round_up_to_even(parsed_inputs["xyz_emmap"].shape[0] * 0.1)
        parsed_inputs["logger"].info("FDR window size is not set. Using a default window size of {} \n".format(fdr_window_size))
    else:
        fdr_window_size = int(parsed_inputs["fdr_w"])

    # Set averaging_filter_size
    averaging_filter_size = parsed_inputs["averaging_filter_size"]    
    
    # Set fdr_filter
    filter_cutoff = choose_fdr_filter_based_on_inputs(parsed_inputs)
    
    fdr_inputs_dict = {
        "emmap_path" : parsed_inputs["xyz_emmap_path"],
        "window_size" : fdr_window_size,
        "fdr" : parsed_inputs["fdr_threshold"],
        "filter_cutoff" : filter_cutoff,
        "averaging_filter_size" : averaging_filter_size,
    }
    
    return fdr_inputs_dict

def prepare_mask_from_inputs(parsed_inputs):
    from locscale.utils.general import get_spherical_mask
    from locscale.preprocessing.headers import run_FDR, check_axis_order
    from locscale.include.emmer.ndimage.map_utils import binarise_map
    from locscale.utils.file_tools import RedirectStdoutToLogger
    if parsed_inputs["mask"] is None:
        # Printing and logging
        print_statement = "A mask path has not been provided. False Discovery Rate control (FDR) based confidence map will be calculated at 1% FDR \n"
        parsed_inputs["logger"].info(print_statement)
        if parsed_inputs["verbose"]:
            tabbed_print.tprint(print_statement)
            
        # Prepare input for FDR
        fdr_inputs_dict = prepare_inputs_for_FDR(parsed_inputs)
        
        # Run FDR with logging
        with RedirectStdoutToLogger(parsed_inputs["logger"], wait_message="Calculating FDR confidence map"):            
            mask_path, mask_path_raw = run_FDR(**fdr_inputs_dict)
            
        xyz_mask_path = check_axis_order(mask_path)
                
        tabbed_print.tprint("Mask generated at {} \n".format(xyz_mask_path))
        parsed_inputs["logger"].info("Mask generated at {} \n".format(xyz_mask_path))
        
        if xyz_mask_path is not None:
            xyz_mask = load_map(xyz_mask_path)[0]
            xyz_mask = binarise_map(xyz_mask, 0.5, return_type='int', threshold_type='gteq')
        else:
            xyz_mask = get_spherical_mask(parsed_inputs["xyz_emmap"].shape)
    else:
        mask_path = parsed_inputs["mask"]
        xyz_mask_path = check_axis_order(mask_path)
        xyz_mask = load_map(xyz_mask_path)[0]
        xyz_mask = binarise_map(xyz_mask,0.5, return_type='int', threshold_type='gteq') # For scaling mask, threshold is always 0.5
        mask_path_raw = mask_path

    return xyz_mask, xyz_mask_path, mask_path_raw

def get_modmap_from_inputs(parsed_inputs):
    from locscale.include.emmer.pdb.pdb_utils import shift_coordinates
    from locscale.utils.file_tools import get_cref_from_inputs
    from locscale.preprocessing.pipeline import get_modmap
    from locscale.include.emmer.ndimage.filter import low_pass_filter
    from locscale.include.emmer.ndimage.map_utils import save_as_mrc, load_map
    
    modality = parsed_inputs["modality"]
    parsed_inputs["logger"].info("Modality: {} \n".format(modality))
    if modality == "map_input_use_directly":
        ## If a model map is provide 
        ## we do not need to use the theoretical profile for scaling
        modmap_path = parsed_inputs["model_map"]
        model_resolution = parsed_inputs["model_resolution"]
        
        ## log the model map path
        parsed_inputs["logger"].info("Model map path: {} \n".format(modmap_path))
        parsed_inputs["logger"].info("Filter resolution: {} \n".format(model_resolution))
        
        if model_resolution is not None:
            if parsed_inputs["verbose"]:
                tabbed_print.tprint("Performing low pass filter on the Model Map \
                    with a cutoff: {} based on user input".format(model_resolution))
            
            pseudo_map_unfiltered_data = load_map(modmap_path)[0]
            pseudo_map_filtered_data = low_pass_filter(im=pseudo_map_unfiltered_data, cutoff=model_resolution, apix=parsed_inputs["apix"])
            
            filename = modmap_path[:-4]+"_filtered.mrc"
            save_as_mrc(map_data=pseudo_map_filtered_data, output_filename=filename, apix=parsed_inputs["apix"])
            
            modmap_path = filename
        
        xyz_modmap_path = check_axis_order(modmap_path)
        xyz_modmap = load_map(xyz_modmap_path)[0]

    else:
        # Collect model map arguments and pass it to get_modmap pipeline
        parsed_inputs["logger"].info("Collecting arguments for model map generation \n")
        pdb_path = parsed_inputs["model_coordinates"]
        
        ## Check if the user has provided the atomic model and set the
        ## use_theoretical_profile to False if yes
        if pdb_path is not None:
            ## If a PDB_path is provided, assume that it is an atomic model hence 
            ## we do not need to use the theoretical profile for scaling
            parsed_inputs["logger"].info("Reading coordinate file: {} \n".format(pdb_path))
            shift_coordinates(in_model_path=pdb_path, trans_matrix=parsed_inputs["shift_vector"],
                                         out_model_path=pdb_path[:-4]+"_shifted.cif")
            pdb_path = pdb_path[:-4]+"_shifted.cif"

    
    
        #############################################################################
        # Stage 4a: Run the get_modmap pipeline                                 #
        #############################################################################
        modmap_path = get_modmap(parsed_inputs)
        xyz_modmap_path = check_axis_order(modmap_path,True)
        xyz_modmap = load_map(xyz_modmap_path)[0]

    return xyz_modmap, xyz_modmap_path     

def set_window_size(window_size_input, apix, verbose):
    from locscale.utils.math_tools import round_up_to_even
    if window_size_input is None:   ## Use default window size of 25 A
        wn = round_up_to_even(25 / apix)
        if verbose:
            tab_print("Using a default window size of {} pixels, corresponding to approximately 25A".format(wn))
    else:
        wn = round_up_to_even(int(window_size_input))
        if verbose:
            tab_print("Provided window size in pixels is {} corresponding to approximately {:.2f} Angstorm".format(wn, wn*apix))
        
    return wn

def pad_input_maps_if_required(parsed_inputs):
    from locscale.utils.general import check_for_window_bleeding, compute_padding_average, pad_or_crop_volume
    parsed_inputs = parsed_inputs.copy()
    wn = parsed_inputs["wn"]
    xyz_emmap = parsed_inputs["xyz_emmap"]
    xyz_mask = parsed_inputs["xyz_mask"]
    xyz_modmap = parsed_inputs["xyz_modmap"]

    window_bleed_and_pad = check_for_window_bleeding(xyz_mask, wn)
    
    parsed_inputs["window_bleed_and_pad"] = window_bleed_and_pad
    ## Collect the padded inputs if required
    if window_bleed_and_pad:
        pad_int_emmap = compute_padding_average(xyz_emmap, xyz_mask)
        pad_int_modmap = compute_padding_average(xyz_modmap, xyz_mask)
        map_shape = [(xyz_emmap.shape[0] + wn), (xyz_emmap.shape[1] + wn), (xyz_emmap.shape[2] + wn)]
        xyz_emmap_new = pad_or_crop_volume(xyz_emmap, map_shape, pad_int_emmap)
        xyz_modmap_new = pad_or_crop_volume(xyz_modmap, map_shape, pad_int_modmap)
        xyz_mask_new = pad_or_crop_volume(xyz_mask, map_shape, 0)
        parsed_inputs["xyz_emmap"] = xyz_emmap_new
        parsed_inputs["xyz_modmap"] = xyz_modmap_new
        parsed_inputs["xyz_mask"] = xyz_mask_new
        parsed_inputs["original_map_shape"] = xyz_emmap.shape
        parsed_inputs["padded_map_shape"] = xyz_emmap_new.shape
    else:
        parsed_inputs["original_map_shape"] = xyz_emmap.shape
        parsed_inputs["padded_map_shape"] = xyz_emmap.shape
    
    return parsed_inputs

def get_scale_factor_arguments(parsed_inputs):
    from locscale.include.emmer.pdb.pdb_tools import find_wilson_cutoff
    from locscale.include.emmer.ndimage.map_tools import compute_radial_profile_simple
    from locscale.include.emmer.ndimage.profile_tools import frequency_array, number_of_segments, estimate_bfactor_through_pwlf

    wilson_cutoff = find_wilson_cutoff(mask_path=parsed_inputs["mask_path_raw"], verbose=False)
    
    resolution_given = parsed_inputs["ref_resolution"] is not None
    poor_resolution = parsed_inputs["ref_resolution"] >= 6 if resolution_given else True
    if (not resolution_given) or poor_resolution:
        high_frequency_cutoff = wilson_cutoff
        nyquist = (round(2*parsed_inputs["apix"]*10)+1)/10
        #fsc_cutoff = fsc_resolution
        bfactor_info = [0,np.array([0,0,0]),np.array([0,0,0])]
        pwlf_fit_quality = 0
    else:
        rp_emmap = compute_radial_profile_simple(parsed_inputs["xyz_emmap"])
        freq = frequency_array(amplitudes=rp_emmap, apix=parsed_inputs["apix"])
        num_segments = number_of_segments(parsed_inputs["ref_resolution"])
        bfactor, amp, (fit,z,slope) = estimate_bfactor_through_pwlf(
            freq=freq, amplitudes=rp_emmap, wilson_cutoff=wilson_cutoff, 
            fsc_cutoff=parsed_inputs["ref_resolution"],num_segments=num_segments, standard_notation=True)
        
        nyquist = (round(2*parsed_inputs["apix"]*10)+1)/10
        #fsc_cutoff = fsc_resolution
        high_frequency_cutoff = 1/np.sqrt(z[-2])
        bfactor_info = [round(bfactor,2), 1/np.sqrt(z).round(2), np.array(slope).round(2)]  ## For information at end
        pwlf_fit_quality = fit.r_squared()
    

    ###############################################################################
    # Stage 6a: Pack into a dictionary
    ###############################################################################

    scale_factor_arguments = {}
    scale_factor_arguments['wilson'] = wilson_cutoff
    scale_factor_arguments['high_freq'] = high_frequency_cutoff
    scale_factor_arguments['fsc_cutoff'] = parsed_inputs["ref_resolution"]
    scale_factor_arguments['nyquist'] = nyquist
    scale_factor_arguments['smooth'] = parsed_inputs["smooth_factor"]
    scale_factor_arguments['boost_secondary_structure'] = parsed_inputs["boost_secondary_structure"]
    scale_factor_arguments['no_reference'] = parsed_inputs["no_reference"]
    scale_factor_arguments['processing_files_folder'] = parsed_inputs["processing_files_folder"]
    scale_factor_arguments['bfactor_info'] = bfactor_info
    scale_factor_arguments['pwlf_fit_quality'] = pwlf_fit_quality
    if parsed_inputs["no_reference"]:
        scale_factor_arguments['set_local_bfactor'] = parsed_inputs["set_local_bfactor"]
    
    return scale_factor_arguments