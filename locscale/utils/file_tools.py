
## FILE HANDLING FUNCTIONS

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

    
def copy_file_to_folder(full_path_to_file, new_folder):
    import shutil
    import os
    
    source = full_path_to_file
    file_name = os.path.basename(source)
    destination = os.path.join(new_folder, file_name)
    shutil.copyfile(source, destination)
    
    return destination

def change_directory(args, folder_name):
    import os    
    from locscale.utils.file_tools import copy_file_to_folder
    
    if folder_name == "processing_files":
        current_directory = os.getcwd()
        new_directory = os.path.join(current_directory, folder_name)
    else:
        new_directory = folder_name
    
    if not os.path.isdir(new_directory):
        os.mkdir(new_directory)
    
    if args.verbose:
        print("Copying files to {}\n".format(new_directory))

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
    
    if args.verbose:
        print("Changing active directory to: {}\n".format(new_directory))
    os.chdir(new_directory)
    
    return args

def generate_filename_from_halfmap_path(in_path):
    ## find filename in the path    
    filename = in_path.split("/")[-1]
    
    ## find EMDB ID in filename
    
    possible_emdb_id = [filename[x:x+4] for x in range(len(filename)-3) if filename[x:x+4].isnumeric()]
    if len(possible_emdb_id) == 1:
        emdb_id = possible_emdb_id[0]
        newfilename = ["EMD_"+emdb_id+"_unfiltered.mrc"]
    else:
        newfilename = ["emdb_map_unfiltered.mrc"]
    
    
    new_path = "/".join(in_path.split("/")[:-1]+newfilename)
    
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
        
        ## TBC

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
    return 
    if args.dev_mode:
        print("Warning: You are in Dev mode. Not checking user input! Results maybe unreliable")
        return 
    
    import mrcfile
    print("ignore profiles default", args.ignore_profiles)
    ## Check input files
    emmap_absent = True
    if args.emmap_path is not None:
        if is_input_path_valid([args.emmap_path]):
            emmap_absent = False
    
    half_maps_absent = True
    if args.half_map1 is not None and args.half_map2 is not None:
        if is_input_path_valid([args.half_map1, args.half_map2]):
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
    
    ## If neither model map or model coordinates are provided, then users cannot use --ignore_profiles and --skip_refine flags
    if model_coordinates_absent and model_map_absent:
        if args.ignore_profiles:
            raise UserWarning("You have not provided a Model Map or Model Coordinates. Thus, pseudo-atomic model will be used for \
                              local sharpening. Please do not raise the --ignore_profiles flag")
        if args.skip_refine:
            raise UserWarning("You have not provided a Model Map or Model Coordinates. Performing REFMAC refinement is essential for \
                              succesful operation of the procedue. Please do not raise the --skip_refine flag")
        
        if args.ref_resolution is None:
            raise UserWarning("You have not provided a Model Map or Model Coordinates. To use REFMAC refinement, resolution target is necessary. \
                              Please provide a target resolution using -res or --ref_resolution")
                            
    
    if model_coordinates_present and model_map_absent:
        if args.skip_refine:
            print("You have asked to skip REFMAC refinement. Atomic bfactors from the input model will be used for simulating Model Map")
        else:
            if args.ref_resolution is None:
                raise UserWarning("You have provided Model Coordinates. By default, the model bfactors will be refined using REFMAC. \
                                  For this, a target resolution is required. Please provide this resolution target using -res or --ref_resolution. \
                                      Instead if you think model bfactors are accurate, then raise the --skip_refine flag to ignore bfactor refinement.")
            

   
    
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
                  




