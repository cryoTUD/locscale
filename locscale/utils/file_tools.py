
## FILE HANDLING FUNCTIONS

from posixpath import isabs


def check_dependencies():
    
    import warnings
    import os
    dependency = {}
    
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
    
        # Check module locscale
        try:
            import locscale
            dependency["locscale"] = True
        except ImportError:
            dependency["locscale"] = False
        
        # Check module gemmi
        try:
            import gemmi
            dependency["gemmi"] = True
        except ImportError:
            dependency["gemmi"] = False
        
        ## Check modules mrcfile, pandas, scipy, numpy, matplotlib, tqdm
        try:
            import mrcfile
            dependency["mrcfile"] = True
        except ImportError:
            dependency["mrcfile"] = False

        # Check module pandas
        try:
            import pandas
            dependency["pandas"] = True
        except ImportError:
            dependency["pandas"] = False
        
        # Check module scipy
        try:
            import scipy
            dependency["scipy"] = True
        except ImportError:
            dependency["scipy"] = False
        
        # Check module numpy
        try:
            import numpy
            dependency["numpy"] = True
        except ImportError:
            dependency["numpy"] = False
        
        # Check module matplotlib
        try:
            import matplotlib
            dependency["matplotlib"] = True
        except ImportError:
            dependency["matplotlib"] = False
        
        ## Check module tqdm
        try:
            import tqdm
            dependency["tqdm"] = True
        except ImportError:
            dependency["tqdm"] = False
        
        ## Check modules tensorflow, keras, tensorflow_addons, pypdb, pyfiglet, emda, proshade
        try:
            import tensorflow
            dependency["tensorflow"] = True
        except:
            dependency["tensorflow"] = False
        
        try:
            import keras
            dependency["keras"] = True
        except:
            dependency["keras"] = False
        
        try:
            import tensorflow_addons
            dependency["tensorflow_addons"] = True
        except:
            dependency["tensorflow_addons"] = False
        
        try:
            import pypdb
            dependency["pypdb"] = True
        except:
            dependency["pypdb"] = False
        
        try:
            import pyfiglet
            dependency["pyfiglet"] = True
        except:
            dependency["pyfiglet"] = False
        
        # try:
        #     import emda
        #     dependency["emda"] = True
        # except:
        #     dependency["emda"] = False
        
        # try:
        #     import proshade
        #     dependency["proshade"] = True
        # except:
        #     dependency["proshade"] = False
        
        ## Check Bio
        try:
            import Bio
            dependency["Bio"] = True
        except:
            dependency["Bio"] = False
        
        ## Check Bio.PDB
        try:
            import Bio.PDB
            dependency["Bio.PDB"] = True
        except:
            dependency["Bio.PDB"] = False
    
    list_of_all_imports = [x for x in dependency.values()]
    if all(list_of_all_imports):
        return True
    else:
        missing_imports = [x for x in dependency.keys() if not dependency[x]]
        return missing_imports

        
def get_locscale_path():
    import locscale
    import os
    return os.path.dirname(locscale.__path__[0])

def get_input_file_directory(args):
    import os
    # Check if halfmap paths are given or emmap path is given
    if args.halfmap_paths is not None:
        halfmap_given = True
    else:
        halfmap_given = False
    
    if halfmap_given:
        # Get the emmap folder from the input paths
        assert os.path.exists(args.halfmap_paths[0]), "Halfmap 1 path does not exist"

        halfmap_1_path_full = os.path.abspath(args.halfmap_paths[0])
        input_folder = os.path.dirname(halfmap_1_path_full)
    else:
        assert os.path.exists(args.emmap_path), "EMmap path does not exist"
        input_folder = os.path.dirname(os.path.abspath(args.emmap_path))

    return input_folder

    
def copy_file_to_folder(full_path_to_file, new_folder):
    import shutil
    import os
    
    source = full_path_to_file
    file_name = os.path.basename(source)
    destination = os.path.join(new_folder, file_name)
    if not os.path.exists(destination):
        shutil.copyfile(source, destination)
    else:
        print(f"File {destination} already exists")
    
    return destination

def change_directory(args, folder_name):
    import os    
    from locscale.utils.file_tools import copy_file_to_folder
    
    # Get the input folder
    input_folder = get_input_file_directory(args)
    
    if folder_name is None:
        new_directory = os.path.join(input_folder, "processing_files")
    else:
        if os.path.isabs(folder_name):
            new_directory = folder_name
        else:
            new_directory = os.path.join(input_folder, folder_name)
    
    if not os.path.isdir(new_directory):
        os.mkdir(new_directory)
    
    assert os.path.isdir(new_directory), "New directory does not exist"
    assert os.path.isabs(new_directory), "New directory is not absolute"
    
    if args.verbose:
        print("Copying files to {}\n".format(new_directory))
    
    # Set the "output_processing_files" argument to the new_directory
    setattr(args, "output_processing_files", new_directory)

    for arg in vars(args):
        value = getattr(args, arg)
        if isinstance(value, str):
            if os.path.exists(value) and arg != "outfile" and arg != "output_processing_files":
                new_location=copy_file_to_folder(value, new_directory)
                setattr(args, arg, new_location)
        if isinstance(value, list):
            if arg == "halfmap_paths":
                halfmap_paths = value
                halfmap1_path = halfmap_paths[0]
                halfmap2_path = halfmap_paths[1]

                new_halfmap1_path = copy_file_to_folder(halfmap1_path, new_directory)
                new_halfmap2_path = copy_file_to_folder(halfmap2_path, new_directory)
                new_halfmap_paths = [new_halfmap1_path,new_halfmap2_path]
                setattr(args, arg, new_halfmap_paths)
    
    return args

def generate_filename_from_halfmap_path(in_path):
    ## find filename in the path    
    import os
    filename = os.path.basename(in_path)

    ## Find all the numbers in the filename
    import re
    numbers = re.findall('[0-9]+',filename)
    
    ## Select only those numbers which have four or five digits
    emdb_numbers = [int(x) for x in numbers if len(str(x)) in [4,5]]

    ## If there is only one number in the filename, then it is the EMDB number
    if len(emdb_numbers) == 1:
        emdb_number = emdb_numbers[0]
        newfilename = "EMD_{}_unsharpened_fullmap.mrc".format(emdb_number)
    else:
        newfilename = "emdb_map_unfiltered.mrc"
    
    new_path = os.path.join(os.path.dirname(in_path), newfilename)
    
    return new_path
    
def get_emmap_path_from_args(args):
    from locscale.utils.file_tools import generate_filename_from_halfmap_path
    from locscale.include.emmer.ndimage.map_tools import add_half_maps
    from locscale.utils.general import shift_map_to_zero_origin
    
    if args.emmap_path is not None:    
        emmap_path = args.emmap_path
        shift_vector=shift_map_to_zero_origin(emmap_path)
    elif args.halfmap_paths is not None:
        print("Adding the two half maps provided to generate a full map \n")
        halfmap_paths = args.halfmap_paths
        assert len(halfmap_paths) == 2, "Please provide two half maps"
        print(halfmap_paths[0])
        print(halfmap_paths[1])

        halfmap1_path = halfmap_paths[0]
        halfmap2_path = halfmap_paths[1]
        new_file_path = generate_filename_from_halfmap_path(halfmap1_path)
        emmap_path = add_half_maps(halfmap1_path, halfmap2_path,new_file_path)
        shift_vector=shift_map_to_zero_origin(halfmap1_path)
    
    return emmap_path, shift_vector
        
def set_modality_based_on_input(args):
    '''
    This function compares all the different ways to get a model map and returns the best option based on the inputs provided
    Different cases available are: 
    1) Complete atomic model provided (or the model map is provided)
    2) Partial atomic model provided
    3) No atomic model provided 
    
    Case 1: Complete atomic model provided
    - if the input is a coordinate file, then refine the model by default and generate a model map
    - if the input is a map file, then use the map file as the model map
    
    Case 2: Partial atomic model provided
    - complete the atomic model using pseudo-atomic model building. Then refine the integrated model and generate model map
    
    Case 3: No atomic model provided
    - Predict the model map using a neural network and return the predicted model map
    - If specified, build a pseudo-atomic model, refine and generate model map   
    '''
    modalities = ["full_atomic_model_refine_and_map", "map_input_use_directly", "partial_atomic_model_build_and_refine", "no_model_predict", "no_model_build_and_refine"]
    
    ## Check input files
    emmap_absent = True
    if args.emmap_path is not None:
        if is_input_path_valid([args.emmap_path]):
            emmap_absent = False
    
    half_maps_absent = True
    if args.halfmap_paths is not None:
        if is_input_path_valid([args.halfmap_paths[0], args.halfmap_paths[1]]):
            half_maps_absent = False
    
    mask_absent = True
    if args.mask is not None:
        if is_input_path_valid([args.mask]):
            mask_absent = False
    
    model_map_absent = True
    if args.model_map is not None:
        if is_input_path_valid([args.model_map]):
            model_map_absent = False
    
    model_coordinates_absent = True
    if args.model_coordinates is not None:
        if is_input_path_valid([args.model_coordinates]):
            model_coordinates_absent = False
    
    emmap_present, half_maps_present = not(emmap_absent), not(half_maps_absent)
    model_map_present, model_coordinates_present = not(model_map_absent), not(model_coordinates_absent)
    
    pseudomodel_not_required = ["map_input_use_directly", "full_model_input_no_refine","full_model_input_refine_and_map", "no_reference","predict_model_map"]
    pseudomodel_required = ["partial_model_input_build_and_refine", "pseudo_model_build_and_refine"]
    if model_map_present:
        modality = "map_input_use_directly"
    elif model_coordinates_present:
        if args.complete_model:
            modality = "partial_model_input_build_and_refine"
        else:
            if args.skip_refine:
                modality = "full_model_input_no_refine"
            else:
                modality = "full_model_input_refine_and_map"
    else:
        if args.no_reference:
            modality = "no_reference"
        else:
            if args.build_using_pseudomodel:
                modality = "pseudo_model_build_and_refine"
            else:
                modality = "predict_model_map"
    
    args.modality = modality
    args.use_theoretical_profile = modality in pseudomodel_required 
    
    print("Running LocScale with modality: {}".format(modality))
    return args
    

def is_input_path_valid(list_of_test_paths):
    '''
    Check if a list of paths are not None and if path points to an actual file

    Parameters
    ----------
    list_of_test_paths : list
        list of paths

    Returns
    -------
    None.

    '''
    import os
    
    for test_path in list_of_test_paths:
        if test_path is None:
            is_test_path_valid = False
            return is_test_path_valid
        if not os.path.exists(test_path):
            is_test_path_valid = False
            return is_test_path_valid
    
    ## If all tests passed then return True
    is_test_path_valid = True
    return is_test_path_valid

def simple_test_model_to_map_fit(args):
    '''
    Test the model to map fit
    '''
    from locscale.include.emmer.pdb.fitting_tools import compute_model_to_map_correlation
    unsharpened_emmap_path, _  = get_emmap_path_from_args(args)
    model_path = args.model_coordinates

    correlation = compute_model_to_map_correlation(emmap_path=unsharpened_emmap_path, pdb_path=model_path)

    return correlation



def check_user_input(args):
    '''
    Check user inputs for errors and conflicts

    Parameters
    ----------
    args : TYPE
        DESCRIPTION.

    Returns
    -------
    None.

    '''
    
    if args.dev_mode:
        print("Warning: You are in Dev mode. Not checking user input! Results maybe unreliable")
        return 
    
    import mrcfile
    
    ## Check input files
    emmap_absent = True
    if args.emmap_path is not None:
        if is_input_path_valid([args.emmap_path]):
            emmap_absent = False
    
    half_maps_absent = True
    if args.halfmap_paths is not None:
        if is_input_path_valid([args.halfmap_paths[0], args.halfmap_paths[1]]):
            half_maps_absent = False
    
    mask_absent = True
    if args.mask is not None:
        if is_input_path_valid([args.mask]):
            mask_absent = False
    
    model_map_absent = True
    if args.model_map is not None:
        if is_input_path_valid([args.model_map]):
            model_map_absent = False
    
    model_coordinates_absent = True
    if args.model_coordinates is not None:
        if is_input_path_valid([args.model_coordinates]):
            model_coordinates_absent = False
    
    ## Rename variables
    emmap_present, half_maps_present = not(emmap_absent), not(half_maps_absent)
    model_map_present, model_coordinates_present = not(model_map_absent), not(model_coordinates_absent)
    ## Sanity checks
    
    ## If emmap is absent or half maps are absent, raise Exceptions
    
    if emmap_absent and half_maps_absent:
        raise UserWarning("Please input either an unsharpened map or two half maps")
          
    
    if model_coordinates_present and model_map_present:
        raise UserWarning("Please provide either a model map or a model coordinates. Not both")
    
    ## If neither model map or model coordinates are provided, then users cannot use --skip_refine flags
    if model_coordinates_absent and model_map_absent:
        warn_against_skip_refine(args, tolerate=False)

        if args.build_using_pseudomodel:
            check_and_warn_about_ref_resolution(args)

                            
    if model_coordinates_present:
        # Check the model to map fit
        correlation = simple_test_model_to_map_fit(args)
        correlation_threshold = 0.3
        if correlation < correlation_threshold:
            print(f"Warning: The model to map correlation is {correlation:.2f}. This is below the threshold of {correlation_threshold}")
        

    if model_coordinates_present and model_map_absent:
        warn_against_skip_refine(args, tolerate=True)
        
    if model_coordinates_present and args.complete_model:
        warn_against_skip_refine(args, tolerate=False)
        check_and_warn_about_ref_resolution(args)
        
    ## Check for window size < 10 A
    if args.window_size is not None:
        window_size_pixels = int(args.window_size)
        if window_size_pixels%2 > 0:
            print("You have input an odd window size. For best performance, an even numbered window size is required. Adding 1 to the provided window size ")
        if args.apix is not None:
            apix = float(args.apix)
        else:
            if args.emmap_path is not None:
                apix = mrcfile.open(args.emmap_path).voxel_size.x
            elif args.half_map1 is not None:
                apix = mrcfile.open(args.half_map1).voxel_size.x
        
        window_size_ang = window_size_pixels * apix
        
        
        if window_size_ang < 10:
            print("Warning: Provided window size of {} is too small for pixel size of {}. \
                  Default window size is generally 25 A. Think of increasing the window size".format(window_size_pixels, apix))
    
    ## Check if the user added a no_reference flag

    if args.no_reference:
        from textwrap import fill
        disclaimer = " Warning: You have asked to not use a reference to perform sharpening. This is not recommended for the following reason. \n\
            With no reference, it is possible to perform local sharpening by scaling the local b-factor to a constant value. \n \
            The constant b-factor is set to 20 by default. When using with reconstruction where the local resolution has a spatial variation \
            then it is not optimal to set the local b-factor to a constant value as it likely boost noise present in high resolution regions. " 
        
        print(fill(disclaimer, width=80))
        ## Pause 
        import time
        time.sleep(2)
        
    ## Find the modalities to generate reference map based on user inputs 

def warn_against_skip_refine(args, tolerate):
    if args.skip_refine:
        if not tolerate:
            raise UserWarning("You have asked to skip REFMAC refinement. \
                                However, you have asked to complete a partially built model. This requires a refined pseudo-atomic model. \
                                Please do not raise the --skip_refine flag")

        if tolerate: 
            print("Warning: You have asked to skip REFMAC refinement. \
                    Please make sure that the atomic ADPs are refined. LocScale performance maybe severely affected if the ADPs are not refined")
                
def check_and_warn_about_ref_resolution(args):
    if args.ref_resolution is None:
        raise UserWarning("Please provide a reference resolution using the --ref_resolution flag")
                

def get_cref_from_inputs(parsed_inputs):
    from locscale.include.emmer.ndimage.filter import get_cosine_mask
    from locscale.include.emmer.ndimage.fsc_util import get_fsc_filter
    from locscale.include.emmer.ndimage.map_utils import load_map

    softmask = get_cosine_mask(parsed_inputs["xyz_mask"], 5)
    halfmap_1, apix = load_map(parsed_inputs["halfmap_paths"][0])
    halfmap_2, apix = load_map(parsed_inputs["halfmap_paths"][1])

    cref = get_fsc_filter(halfmap_1*softmask, halfmap_2*softmask)
    return cref

def get_cref_from_arguments(args, mask):
    '''
    Get the cref value from the arguments
    
    Parameters
    ----------
    args : TYPE
        DESCRIPTION.

    Returns
    -------
    cref : TYPE
        DESCRIPTION.

    '''
    from locscale.include.emmer.ndimage.fsc_util import get_fsc_filter
    from locscale.include.emmer.ndimage.map_utils import load_map
    ## Check if halfmaps present in arguments
    if args.halfmap_paths is not None:
        half_maps_present = True
    else:
        half_maps_present = False
    
    if half_maps_present:
        halfmap_path_1 = args.halfmap_paths[0]
        halfmap_path_2 = args.halfmap_paths[1]

        halfmap_1, apix = load_map(halfmap_path_1)
        halfmap_2, apix = load_map(halfmap_path_2)

        masked_halfmap_1 = halfmap_1 * mask
        masked_halfmap_2 = halfmap_2 * mask
        Cref = get_fsc_filter(masked_halfmap_1, masked_halfmap_2)
        
    else:
        Cref = None
    
    return Cref

def get_fsc_curve_from_arguments(args):
    '''
    Get the fsc curve from the arguments
    
    Parameters
    ----------
    args : TYPE
        DESCRIPTION.

    Returns
    -------
    fsc_curve : TYPE
        DESCRIPTION.

    '''
    from locscale.include.emmer.ndimage.fsc_util import calculate_fsc_maps
    ## Check if halfmaps present in arguments
    if args.halfmap_paths is not None:
        half_maps_present = True
    else:
        half_maps_present = False
    
    if half_maps_present:
        halfmap_path_1 = args.halfmap_paths[0]
        halfmap_path_2 = args.halfmap_paths[1]

        fsc_curve = calculate_fsc_maps(halfmap_path_1, halfmap_path_2)
    else:
        fsc_curve = None
    
    return fsc_curve

        
    
                  




