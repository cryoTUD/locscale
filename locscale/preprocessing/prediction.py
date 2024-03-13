
import os
import numpy as np
from locscale.utils.plot_tools import tab_print
from locscale.utils.file_tools import RedirectStdoutToLogger, pretty_print_dictionary
from locscale.emmernet.run_emmernet import run_emmernet
from locscale.include.emmer.ndimage.map_utils import load_map, save_as_mrc

tabbed_print = tab_print(2)
tprint = tabbed_print.tprint

def predict_model_map_from_input_map(parsed_inputs):
    from locscale.include.emmer.ndimage.map_utils import load_map, save_as_mrc
    from locscale.emmernet.utils import check_emmernet_dependencies, check_and_download_emmernet_model
    emmernet_model_folder = check_and_download_emmernet_model(verbose=True)
    
    # set inputs from parsed_inputs
    cube_size = 32 
    symmetry = parsed_inputs["symmetry"]
    emmap_path = parsed_inputs["xyz_emmap_path"]
    xyz_mask_path = parsed_inputs["mask_path_raw"]
    if parsed_inputs["use_low_context_model"]:
        trained_model = "emmernet_low_context"
    else:
        trained_model = "emmernet_high_context"
    stride = 16 
    batch_size = parsed_inputs["batch_size"]
    gpu_ids = parsed_inputs["gpu_ids"]
    verbose = parsed_inputs["verbose"]
    target_map_path = None
    model_path = parsed_inputs["model_path"]
    monte_carlo = False
    monte_carlo_iterations = 1
    physics_based = False
    
    input_dictionary = {}
    input_dictionary["cube_size"] = cube_size
    input_dictionary["emmap_path"] = emmap_path
    input_dictionary["xyz_mask_path"] = xyz_mask_path
    input_dictionary["trained_model"] = trained_model
    input_dictionary["stride"] = stride
    input_dictionary["batch_size"] = batch_size
    input_dictionary["gpu_ids"] = gpu_ids
    input_dictionary["verbose"] = verbose
    input_dictionary["emmernet_model_folder"] = emmernet_model_folder
    input_dictionary["target_map_path"] = target_map_path
    input_dictionary["model_path"] = model_path
    input_dictionary["monte_carlo"] = monte_carlo
    input_dictionary["monte_carlo_iterations"] = monte_carlo_iterations
    input_dictionary["physics_based"] = physics_based
    input_dictionary["logger"] = parsed_inputs["logger"]
    input_dictionary["symmetry"] = parsed_inputs["symmetry"]
    input_dictionary["output_processing_files"] = parsed_inputs["output_processing_files"]
    if gpu_ids is None:
        cuda_visible_devices_string = ""
    else:
        cuda_visible_devices_string = ",".join([str(gpu_id) for gpu_id in gpu_ids])
    
    os.environ["CUDA_VISIBLE_DEVICES"] = cuda_visible_devices_string
    if verbose:
        print("\tCUDA_VISIBLE_DEVICES: {}".format(os.environ["CUDA_VISIBLE_DEVICES"]))

    parsed_inputs["logger"].info("Input dictionary: \n{}".format(pretty_print_dictionary(input_dictionary)))
    
    # run emmernet
    emmap, apix = load_map(emmap_path)
    input_dictionary["apix"] = apix
    emmernet_output = run_emmernet(input_dictionary)
    model_map_predicted = emmernet_output["output_predicted_map_mean"]
    emmap_extension = os.path.splitext(emmap_path)[1]
    model_map_path_filename = emmap_path.replace(emmap_extension, "_model_map_predicted.mrc")
    save_as_mrc(model_map_predicted, model_map_path_filename, apix)
    # if symmetry != "C1":
    #     from locscale.include.symmetry_emda.symmetrize_map import symmetrize_map_emda
    #     from locscale.include.emmer.ndimage.map_utils import save_as_mrc
        
    #     if verbose:
    #         print_statement = "b) Applying symmetry: {}".format(symmetry)
    #         print(print_statement)
    #         parsed_inputs['logger'].info(print_statement)
        
    #     with RedirectStdoutToLogger(parsed_inputs['logger'], wait_message="Applying symmetry"):
    #         sym = symmetrize_map_emda(emmap_path=model_map_path_filename,pg=symmetry)
    #         symmetrised_modmap = model_map_path_filename[:-4]+"_{}_symmetry.mrc".format(symmetry)
    #         save_as_mrc(map_data=sym, output_filename=symmetrised_modmap, apix=apix, origin=0, verbose=True)
    #         predicted_modmap = symmetrised_modmap
    # else:
    #     predicted_modmap = model_map_path_filename
    #     if verbose:
    #         print_statement = "b) No symmetry applied"
    #         print(print_statement)
    #         parsed_inputs['logger'].info(print_statement)
    
    print_statement = "Predicted model map with shape {} saved to: {}".format(model_map_predicted.shape, model_map_path_filename)
    parsed_inputs["logger"].info(print_statement)
    tprint(print_statement)
    
    
    return model_map_path_filename
    