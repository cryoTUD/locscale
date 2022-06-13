## Script to run EMmerNet on an input map
## import the necessary packages from locscale.include.emmer

from locscale.include.emmer.ndimage.map_utils import resample_map, load_map
from locscale.emmernet.emmernet_functions import standardize_map, minmax_normalize_map, get_cubes, assemble_cubes
import numpy as np
import os

os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3'  


def run_emmernet(input_dictionary):
    ## Ignore DeprecationWarning
    import warnings
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", category=DeprecationWarning)
        import tensorflow as tf
        from tensorflow.keras.models import load_model
    
    from locscale.utils.general import merge_sequence_of_sequences
    EMMERNET_CUBE_SIZE=32
    
    ## Get the input map path
    emmap_path = input_dictionary["emmap_path"]
    emmernet_type = input_dictionary["emmernet_type"]
    stride = input_dictionary["stride"]
    batch_size = input_dictionary["batch_size"]
    gpu_ids = input_dictionary["gpu_ids"]
    verbose = input_dictionary["verbose"]
    emmernet_model_folder = input_dictionary["emmernet_model_folder"]

    emmap, apix = load_map(emmap_path)
    if verbose:
        print("Emmap loaded from: {}".format(emmap_path))
        print("Emmap shape: {}".format(emmap.shape))
        print("Pixelsize read as: {:.2f}".format(apix))

        print("1) Pre-processing commencing...")
    
    ## Preprocess

    emmap_preprocessed = preprocess_map(emmap, apix)
    if verbose:
        print("\tPreprocessing complete")
        print("\tPre-processed map shape: {}".format(emmap_preprocessed.shape))
        print("2) Prediction commencing...")

    cubes, cubecenters = get_cubes(emmap_preprocessed, cube_size=EMMERNET_CUBE_SIZE, step_size=stride)
    if verbose:
        print("\tCubes extracted")
        print("\tNumber of cubes: {}".format(len(cubes)))
    ## Load the model

    emmernet_model = load_emmernet_model(emmernet_type, emmernet_model_folder)
    if verbose:
        print("\tEMmerNet model loaded: {}".format(emmernet_type))

    ## Run EMmerNet using GPUs
    
    # prepare GPU id list
    if gpu_ids is None:
        print("No GPU id specified, running on CPU")
        print("If you want to use GPUs, please specify the GPU id(s) using the --gpu_ids flag")
        print("This may take a while...")
        mirrored_strategy = tf.distribute.MirroredStrategy()
    else:
        gpu_id_list = ["/gpu:"+str(gpu_id) for gpu_id in gpu_ids]
        if verbose:
            print("\tGPU ids: {}".format(gpu_id_list))
        mirrored_strategy = tf.distribute.MirroredStrategy(devices=gpu_id_list)
    
    predicted_cubes = run_emmernet_batch(cubes, emmernet_model, mirrored_strategy, batch_size=batch_size)

    if verbose:
        print("\tEMmerNet prediction complete")
        print("\tNumber of predicted cubes: {}".format(len(predicted_cubes)))
    ## Merge the predicted cubes sequence
    
    predicted_map = assemble_cubes(predicted_cubes, cubecenters, emmap_preprocessed.shape[0],EMMERNET_CUBE_SIZE)
    if verbose:
        print("\tPredicted map assembled")
        print("\tPredicted map shape: {}".format(predicted_map.shape))
        print("3) Post-processing commencing...")
    
    
    ## Postprocess

    predicted_map_postprocessed = postprocess_map(predicted_map, apix)
    if verbose:
        print("\tPost-processing complete")
        print("\tPost-processed map shape: {}".format(predicted_map_postprocessed.shape))

    #return predicted_map_postprocessed

    emmernet_output_dictionary = {"output":predicted_map_postprocessed}

    return emmernet_output_dictionary


def load_emmernet_model(emmernet_type, emmernet_model_folder=None):
    import os
    ## Ignore DeprecationWarning
    import warnings

    if emmernet_model_folder is None:
        import locscale
        emmernet_model_folder = os.path.join(os.path.dirname(locscale.__file__), "emmernet", "emmernet_models")
    
    assert os.path.exists(emmernet_model_folder), "EMmerNet model folder not found: {}".format(emmernet_model_folder)

    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", category=DeprecationWarning)    
        from tensorflow.keras.models import load_model
        from tensorflow_addons.layers import GroupNormalization

    emmernet_folder_path = emmernet_model_folder
    if emmernet_type == "model_based":
        emmernet_model_path = os.path.join(emmernet_folder_path, "EMmerNet_MB.hdf5")
    elif emmernet_type == "model_free":
        emmernet_model_path = os.path.join(emmernet_folder_path, "EMmerNet_MF.hdf5")
    elif emmernet_type == "ensemble":
        emmernet_model_path = os.path.join(emmernet_folder_path, "EMmerNet_MBMF.hdf5")
    else:
        raise ValueError("Invalid emmernet_type")
    
    emmernet_model = load_model(emmernet_model_path)
    
    return emmernet_model

def run_emmernet_batch(cubes, emmernet_model, mirrored_strategy, batch_size):
    ## Run the model on the cube
    ## Ignore DeprecationWarning
    import warnings
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", category=DeprecationWarning)
        import tensorflow_datasets as tfds
        import atexit
    from tqdm import tqdm
    import os
    

    tfds.disable_progress_bar()
    cube_size = cubes[0].shape[0]
    cubes = np.array(cubes)
    cubes_x = np.expand_dims(cubes, axis=4)
    cubes_predicted = np.empty((0, cube_size, cube_size, cube_size, 1))
    with mirrored_strategy.scope():
        for i in tqdm(np.arange(0,len(cubes),batch_size),desc="Running EMmerNet"):
            cubes_batch_X = np.empty((batch_size, cube_size, cube_size, cube_size, 1))
            cubes_batch_X = cubes_x[i:i+batch_size,:,:,:,:]

            cubes_batch_predicted = emmernet_model.predict(x=cubes_batch_X, batch_size=batch_size, verbose=0)

            cubes_predicted = np.append(cubes_predicted, cubes_batch_predicted, axis=0)
        
        
    # close the mirrored strategy's multiprocessing ThreadPool explicitly
    atexit.register(mirrored_strategy._extended._collective_ops._pool.close)

    # squeeze cubes to 3 dimensions
    cubes_predicted = np.squeeze(cubes_predicted, axis=-1)
    return cubes_predicted


## Preprocess the map
def preprocess_map(emmap, apix):
    ## Resample the map to 1A per pixel
    emmap_resampled = resample_map(emmap, apix=apix,apix_new=1)

    ## standardize the map
    emmap_standardized = standardize_map(emmap_resampled)

    return emmap_standardized

def postprocess_map(predicted_map, apix):
    ## Resample the map to the original pixel size
    predicted_map_resampled = resample_map(predicted_map, apix=1,apix_new=apix)

    ## MinMax normalize the map
    predicted_map_normalized = minmax_normalize_map(predicted_map_resampled)

    return predicted_map_normalized



    