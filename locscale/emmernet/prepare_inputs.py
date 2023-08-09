import mrcfile
import os
import numpy as np
def prepare_inputs(args):
    import os
    from locscale.emmernet.utils import check_emmernet_dependencies, check_and_download_emmernet_model
    from locscale.utils.file_tools import get_emmap_path_from_args, pretty_print_dictionary
    from locscale.include.emmer.ndimage.map_utils import load_map
    from locscale.preprocessing.headers import check_axis_order
    from locscale.utils.prepare_inputs import prepare_mask_from_inputs
    print("."*80)

    check_emmernet_dependencies(verbose=True)
    emmernet_model_folder = check_and_download_emmernet_model(verbose=True)
    
    parsed_inputs = vars(args)
    ###########################################################################
    # Stage 1: Extract the map
    ###########################################################################
    parsed_inputs["unsharpened_emmap_path"], parsed_inputs["shift_vector"]  = get_emmap_path_from_args(args)
    # Check axis orders of the maps
    parsed_inputs["xyz_emmap_path"] = check_axis_order(parsed_inputs["unsharpened_emmap_path"])
    parsed_inputs["xyz_emmap"], apix_from_file = load_map(parsed_inputs["xyz_emmap_path"])
    parsed_inputs["apix"] = apix_from_file
    parsed_inputs["emmap_path"] = parsed_inputs["xyz_emmap_path"]
    parsed_inputs["emmap_folder"] = os.path.dirname(parsed_inputs["emmap_path"])
    # Monte Carlo 
    parsed_inputs["monte_carlo"] = not(parsed_inputs["no_monte_carlo"])
    
    # GPU IDs
    gpu_ids = parsed_inputs["gpu_ids"]
    
    # Verbose
    verbose = parsed_inputs["verbose"]
    
    # Model to use
    if parsed_inputs["use_low_context_model"]:
        trained_model = "emmernet_low_context"
    else:
        trained_model = "emmernet_high_context"
    
    parsed_inputs["trained_model"] = trained_model
    
    if gpu_ids is None:
        cuda_visible_devices_string = ""
    else:
        cuda_visible_devices_string = ",".join([str(gpu_id) for gpu_id in gpu_ids])
    
    os.environ["CUDA_VISIBLE_DEVICES"] = cuda_visible_devices_string
    if verbose:
        print("\tCUDA_VISIBLE_DEVICES set to {}".format(os.environ["CUDA_VISIBLE_DEVICES"]))
    
    ###########################################################################
    # Stage 2: Prepare the mask
    ###########################################################################
    
    if parsed_inputs["verbose"]:
        print("."*80)
        print("Preparing mask \n")
    
    parsed_inputs["xyz_mask"], parsed_inputs["xyz_mask_path"], parsed_inputs["mask_path_raw"] = prepare_mask_from_inputs(parsed_inputs)

    parsed_inputs["emmernet_model_folder"] = emmernet_model_folder
    
    if parsed_inputs["verbose"]:
        print("Inputs parsed successfully")
    print("."*80)
    ###########################################################################
    # Print the inputs
    ###########################################################################
    for key, val in parsed_inputs.items():
        if isinstance(val, str):
            print("{}:\t {}".format(key, parsed_inputs[key]))
        elif isinstance(val, bool):
            print("{}:\t {}".format(key, parsed_inputs[key]))
        elif isinstance(val, (list, tuple)):
            print("{}:\t {} (length)".format(key, len(parsed_inputs[key])))
        elif isinstance(val, np.ndarray):
            print("{}:\t {} (shape)".format(key, parsed_inputs[key].shape))
        elif val is not None:
            print("{}:\t {}".format(key, str(parsed_inputs[key])))
        else:
            print("{}:\t --".format(key))
    print("."*80)
    return parsed_inputs







