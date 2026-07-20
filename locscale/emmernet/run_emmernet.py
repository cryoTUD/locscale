#
# Delft University of Technology (TU Delft) hereby disclaims all copyright interest in the program 'LocScale'
# written by the Author(s).
# Copyright (C) 2021 Alok Bharadwaj and Arjen J. Jakobi
# This software may be modified and distributed under the terms of the BSD license. 
# You should have received a copy of the BSD 3-clause license along with this program (see LICENSE file file for details).
# If not see https://opensource.org/license/bsd-3-clause/.
#

## Script to run EMmerNet on an input map
## import the necessary packages from locscale.include.emmer
import os 
import numpy as np 
from scipy.stats import norm 

from locscale.include.emmer.ndimage.map_utils import resample_map, load_map
from locscale.utils.file_tools import RedirectStdoutToLogger
from locscale.emmernet.emmernet_functions import standardize_map, get_cubes, assemble_cubes, replace_cubes_in_dictionary,\
                                                    load_smoothened_mask, show_signal_cubes

from locscale.emmernet.utils import symmetrise_if_needed             

from .emmernet_torch import create_dataloader

def run_emmernet(input_dictionary):
    input_dictionary["logger"].info("1) Preprocessing the data...")
    input_dictionary = start_preprocessing_data(input_dictionary)
    
    input_dictionary["logger"].info("2) Preparing inputs for the network...")
    input_dictionary = prepare_inputs_for_network(input_dictionary)
    
    input_dictionary["logger"].info("3) Predicting the cubes...")
    output_dictionary = predict_cubes_and_assemble(input_dictionary)
    
    output_dictionary = symmetrise_if_needed(input_dictionary=input_dictionary, output_dictionary=output_dictionary)
    
    return output_dictionary

def start_preprocessing_data(input_dictionary):
    emmap_path = input_dictionary["emmap_path"]
    mask_path = input_dictionary["xyz_mask_path"]
    verbose = input_dictionary["verbose"]
    processing_files_folder = os.path.dirname(emmap_path)
    
    with RedirectStdoutToLogger(input_dictionary["logger"], wait_message="Loading input"):
        emmap, apix = load_map(emmap_path, verbose=True)
    
    with RedirectStdoutToLogger(input_dictionary["logger"], wait_message="Loading mask"):
        mask, _ = load_smoothened_mask(mask_path, verbose=True)
    
    input_map_shape = emmap.shape
    
    emmap_preprocessed = preprocess_map(emmap, apix)
    mask_preprocessed = preprocess_map(mask, apix, standardize=False)
    
    if verbose:
        print("\tPreprocessing complete")
        print_statement = "\tPre-processed map shape: {}".format(emmap_preprocessed.shape)
        print(print_statement)
        input_dictionary["logger"].info(print_statement)        

    input_dictionary["emmap_preprocessed"] = emmap_preprocessed
    input_dictionary["mask_preprocessed"] = mask_preprocessed
    input_dictionary["input_map_shape"] = input_map_shape
    input_dictionary["preprocessed_map_shape"] = emmap_preprocessed.shape
    input_dictionary["apix_raw"] = apix
    input_dictionary["processing_files_folder"] = processing_files_folder
    
    return input_dictionary
    
    
def prepare_inputs_for_network(input_dictionary):
    emmap_preprocessed = input_dictionary["emmap_preprocessed"]
    cube_size = input_dictionary["cube_size"]
    stride = input_dictionary["stride"]
    mask_preprocessed = input_dictionary["mask_preprocessed"]
    processing_files_folder = input_dictionary["processing_files_folder"]
    verbose = input_dictionary["verbose"]
    
    cubes_dictionary, cubes_array, signal_cubes = get_cubes(emmap_preprocessed, cube_size=cube_size, step_size=stride, mask=mask_preprocessed)
    cubes_center_save_path = os.path.join(processing_files_folder, "signal_cubes_resampled.mrc")
    show_signal_cubes(signal_cubes, emmap_preprocessed.shape, \
            save_path=cubes_center_save_path, apix=input_dictionary["apix_raw"], input_shape=input_dictionary["input_map_shape"])

    input_dictionary["cubes_dictionary"] = cubes_dictionary
    input_dictionary["cubes_array"] = cubes_array
        
    if verbose:
        print("\tCubes extracted")
        print_statement_cubes = f"\tNumber of cubes: {len(cubes_dictionary)} of which {len(signal_cubes)} are signal cubes"
        print(print_statement_cubes)
        input_dictionary["logger"].info(print_statement_cubes)
        print("\tCheck the centers of the cubes in the file: {}".format(cubes_center_save_path))
        input_dictionary["logger"].info("\tCheck the centers of the cubes in the file: {}".format(cubes_center_save_path))
    return input_dictionary    

def predict_cubes_and_assemble(input_dictionary):
    import os 
    verbose = input_dictionary["verbose"]
    processing_files_folder = input_dictionary["output_processing_files"]
    
    gpu_ids = input_dictionary["gpu_ids"]

    
    if verbose:
        print("\tCUDA_VISIBLE_DEVICES set to: {}".format(os.environ["CUDA_VISIBLE_DEVICES"]))
        input_dictionary["logger"].info("\tCUDA_VISIBLE_DEVICES set to: {}".format(os.environ["CUDA_VISIBLE_DEVICES"]))

    emmernet_model, device = load_emmernet_model(model_type=input_dictionary["trained_model"], verbose=verbose)

    # use nn.DataParallel if multiple GPUs are available
    if gpu_ids is not None and len(gpu_ids) > 1:
        import torch
        emmernet_model = torch.nn.DataParallel(emmernet_model, device_ids=[i for i in range(len(gpu_ids))])
        if verbose:
            print("\tUsing DataParallel on GPUs: {}".format(gpu_ids))
            input_dictionary["logger"].info("\tUsing DataParallel on GPUs: {}".format(gpu_ids))


    input_dictionary["logger"].info("Prediction start")
    input_dictionary = run_emmernet_batch(input_dictionary, emmernet_model, device)

    if input_dictionary["monte_carlo"]:
        input_dictionary["logger"].info("Assembling the cubes in the right place...")
        predicted_map_mean = assemble_cubes_in_right_place(input_dictionary, input_dictionary["cubes_predicted_mean"])
        predicted_map_var = assemble_cubes_in_right_place(input_dictionary, input_dictionary["cubes_predicted_var"])
    else: 
        predicted_map_mean = assemble_cubes_in_right_place(input_dictionary, input_dictionary["cubes_predicted_mean"])
        predicted_map_var = None

    emmernet_output_dictionary = {
        "output_predicted_map_mean":predicted_map_mean, 
        "output_predicted_map_var":predicted_map_var,
        "output_processing_files" : processing_files_folder,
    }
    
    return emmernet_output_dictionary

def assemble_cubes_in_right_place(input_dictionary, predicted_cubes):
        predicted_cubes_dictionary = replace_cubes_in_dictionary(predicted_cubes, input_dictionary["cubes_dictionary"])
        predicted_map_potential = assemble_cubes(predicted_cubes_dictionary,input_dictionary["preprocessed_map_shape"],average=True)
        predicted_map_postprocessed = postprocess_map(predicted_map_potential, input_dictionary["apix_raw"], output_shape=input_dictionary["input_map_shape"])
        
        return predicted_map_postprocessed
   

def get_device():
    import torch
    if torch.cuda.is_available():
        device = torch.device("cuda")
        print("Using GPU: {}".format(torch.cuda.get_device_name(device)))
    elif torch.backends.mps.is_available():
        device = torch.device("mps")
        print("Using Apple Metal Performance Shaders (MPS) backend")
    else:
        device = torch.device("cpu")
        print("Using CPU")
    return device
    
def load_emmernet_model(model_type, verbose=False):
    import os 
    import locscale 
    import torch 
    from .emmernet_torch import EMmerNet

    assert type(model_type) == str, "model_type must be a string"

    default_emmernet_model_folder = os.path.join(os.path.dirname(locscale.__file__), "emmernet", "emmernet_models")
    if model_type == "emmernet_high_context":
        emmernet_model_path = os.path.join(default_emmernet_model_folder, "emmernet", "emmernet_highcontext.pt")
    elif model_type == "emmernet_low_context":
        emmernet_model_path = os.path.join(default_emmernet_model_folder, "emmernet", "emmernet_lowcontext.pt")
    elif os.path.exists(model_type) and model_type.endswith(".pt"):
        emmernet_model_path = model_type
    else:
        raise ValueError("Invalid model_type. Must be 'emmernet_high_context', 'emmernet_low_context' or a valid path to a .pt file.")
    
    device = get_device()
    
    assert os.path.exists(emmernet_model_path), "EMmerNet model file not found at {}".format(emmernet_model_path)

    
    emmernet_model = EMmerNet()

    # load the model weights
    model_weights = torch.load(emmernet_model_path, map_location=device)
    
    # load the model weights into the model
    emmernet_model.load_state_dict(model_weights)

    # send the model to the device
    emmernet_model.to(device)

    # set the model to evaluation mode
    emmernet_model.eval()        
    if verbose:
        print("\tEMmerNet model loaded from: {}".format(emmernet_model_path))
    

    return emmernet_model, device 

    
def run_emmernet_batch(input_dictionary, emmernet_model, device):
    # collect inputs from input dictionary
    monte_carlo = input_dictionary["monte_carlo"]
    monte_carlo_iterations = input_dictionary["monte_carlo_iterations"]
    batch_size = input_dictionary["batch_size"]
    cubes = input_dictionary["cubes_array"]
    cuda_visible_devices_string = input_dictionary["cuda_visible_devices_string"]
    print("Running EMmerNet on {} cubes".format(len(cubes)))
    input_dictionary["logger"].info("Running EMmerNet on {} cubes".format(len(cubes)))
    input_dictionary["logger"].info("Device: {}".format(device))
    input_dictionary["logger"].info("CUDA_VISIBLE_DEVICES: {}".format(cuda_visible_devices_string))

    dataloader = create_dataloader(cubes, batch_size=batch_size, shuffle=False)

    if monte_carlo:
        cubes_predicted_mean, cubes_predicted_var = run_emmernet_pytorch_montecarlo(
            dataloader=dataloader, 
            emmernet_model=emmernet_model,
            device=device,
            num_samples=monte_carlo_iterations
        )

        cubes_predicted_mean = np.squeeze(cubes_predicted_mean, axis=1)
        cubes_predicted_var = np.squeeze(cubes_predicted_var, axis=1)
    else:
        cubes_predicted_mean = run_emmernet_pytorch(
            dataloader=dataloader, 
            emmernet_model=emmernet_model,
            device=device
        )
        cubes_predicted_var = None
    
            
    input_dictionary["cubes_predicted_mean"] = cubes_predicted_mean
    input_dictionary["cubes_predicted_var"] = cubes_predicted_var
    
    return input_dictionary


def run_emmernet_pytorch(dataloader, emmernet_model, device):
    import torch
    from tqdm import tqdm
    import numpy as np

    emmernet_model.eval()
    cubes_predicted = []

    with torch.no_grad():
        for batch in tqdm(dataloader, desc="Running EMmerNet"):
            cubes_batch_X = batch[0].to(device)
            cubes_batch_predicted = emmernet_model(cubes_batch_X)
            cubes_predicted.append(cubes_batch_predicted.cpu().numpy())

    cubes_predicted = np.concatenate(cubes_predicted, axis=0)
    return cubes_predicted

def run_emmernet_pytorch_montecarlo(dataloader, emmernet_model, device, num_samples=15):
    import torch
    import numpy as np
    from tqdm import tqdm

    emmernet_model.train()          # dropout ON; no BatchNorm in EMmerNet, so this
                                    # enables MC dropout and changes nothing else
    means, variances = [], []
    with torch.no_grad():
        for batch in tqdm(dataloader, desc="Running EMmerNet MC"):
            x = batch[0].to(device)
            samples = torch.stack([emmernet_model(x) for _ in range(num_samples)], dim=0)
            var, mean = torch.var_mean(samples, dim=0, unbiased=False)
            means.append(mean.cpu().numpy())
            variances.append(var.cpu().numpy())

    return np.concatenate(means, axis=0), np.concatenate(variances, axis=0)     

## Preprocess the map
def preprocess_map(emmap, apix, standardize=True):
    ## Resample the map to 1A per pixel
    emmap_resampled = resample_map(emmap, apix=apix,apix_new=1)
    ## standardize the map
    if standardize:
        emmap_standardized = standardize_map(emmap_resampled)
        return emmap_standardized
    else:
        return emmap_resampled

def postprocess_map(predicted_map, apix, output_shape):
    ## Resample the map to the original pixel size
    predicted_map_resampled = resample_map(predicted_map, apix=1,apix_new=apix, assert_shape=output_shape)
    return predicted_map_resampled

            
