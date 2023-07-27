
import os
import numpy as np
from locscale.utils.plot_tools import tab_print
from locscale.emmernet.run_emmernet import run_emmernet
from locscale.include.emmer.ndimage.map_utils import load_map, save_as_mrc

tabbed_print = tab_print(2)
tprint = tabbed_print.tprint

def predict_model_map_from_input_map(parsed_inputs):
    from locscale.emmernet.utils import check_emmernet_dependencies, check_and_download_emmernet_model
    tprint("Predicting model map from input map")
    emmernet_model_folder = check_and_download_emmernet_model(verbose=True)
    
    # set inputs from parsed_inputs
    cube_size = 32 
    emmap_path = parsed_inputs["xyz_emmap_path"]
    xyz_mask_path = parsed_inputs["mask_path_raw"]
    if parsed_inputs["prefer_low_context_model"]:
        trained_model = "atomic_model_map_target"
    else:
        trained_model = "hybrid_model_map_target"
    stride = 16 
    batch_size = parsed_inputs["batch_size"]
    gpu_ids = parsed_inputs["gpu_ids"]
    verbose = parsed_inputs["verbose"]
    target_map_path = None
    model_path = parsed_inputs["model_path"]
    monte_carlo = False
    monte_carlo_iterations = 1
    
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

    # run emmernet
    emmap, apix = load_map(emmap_path)
    
    emmernet_output = run_emmernet(input_dictionary)
    model_map_predicted = emmernet_output["output_predicted_map_mean"]
    tprint("Predicted model map shape: {}".format(model_map_predicted.shape))
    emmap_extension = os.path.splitext(emmap_path)[1]
    model_map_path_filename = emmap_path.replace(emmap_extension, "_model_map_predicted.mrc")
    tprint("Saving model map to: {}".format(model_map_path_filename))
    save_as_mrc(model_map_predicted, model_map_path_filename, apix)
    
    return model_map_path_filename
    