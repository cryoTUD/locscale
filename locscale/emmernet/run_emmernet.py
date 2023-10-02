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

from locscale.include.emmer.ndimage.map_utils import resample_map, load_map
from locscale.utils.file_tools import RedirectStdoutToLogger
from locscale.emmernet.emmernet_functions import standardize_map, get_cubes, assemble_cubes, replace_cubes_in_dictionary,\
                                                    load_smoothened_mask, show_signal_cubes

from locscale.emmernet.utils import symmetrise_if_needed                                                    
import tensorflow as tf
import numpy as np
import os
from tqdm import tqdm
from scipy.stats import norm 

os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3'  
tf.random.set_seed(42)

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
    
    verbose = input_dictionary["verbose"]

    mirrored_strategy, cuda_visible_devices_string = get_strategy(input_dictionary)
    os.environ["CUDA_VISIBLE_DEVICES"] = cuda_visible_devices_string
    input_dictionary["cuda_visible_devices_string"] = cuda_visible_devices_string
    if verbose:
        print("\tCUDA_VISIBLE_DEVICES set to: {}".format(os.environ["CUDA_VISIBLE_DEVICES"]))
        input_dictionary["logger"].info("\tCUDA_VISIBLE_DEVICES set to: {}".format(os.environ["CUDA_VISIBLE_DEVICES"]))

    
    emmernet_model = load_emmernet_model(input_dictionary)
    input_dictionary["logger"].info("Prediction start")
    input_dictionary = run_emmernet_batch(input_dictionary, emmernet_model, mirrored_strategy)

    if input_dictionary["monte_carlo"] or input_dictionary["physics_based"]:
        input_dictionary["logger"].info("Assembling the cubes in the right place...")
        predicted_map_mean = assemble_cubes_in_right_place(input_dictionary, input_dictionary["cubes_predicted_mean"])
        predicted_map_var = assemble_cubes_in_right_place(input_dictionary, input_dictionary["cubes_predicted_var"])
        predicted_map_total = assemble_cubes_in_right_place(input_dictionary, input_dictionary["cubes_predicted_total"])
    else: 
        predicted_map_mean = assemble_cubes_in_right_place(input_dictionary, input_dictionary["cubes_predicted_mean"])
        predicted_map_var = None
        predicted_map_total = None

    emmernet_output_dictionary = {
        "output_predicted_map_mean":predicted_map_mean, 
        "output_predicted_map_var":predicted_map_var,
        "output_predicted_map_total":predicted_map_total,
    }

    
    return emmernet_output_dictionary

def assemble_cubes_in_right_place(input_dictionary, predicted_cubes):
        predicted_cubes_dictionary = replace_cubes_in_dictionary(predicted_cubes, input_dictionary["cubes_dictionary"])
        predicted_map_potential = assemble_cubes(predicted_cubes_dictionary,input_dictionary["preprocessed_map_shape"],average=True)
        predicted_map_postprocessed = postprocess_map(predicted_map_potential, input_dictionary["apix_raw"], output_shape=input_dictionary["input_map_shape"])
        
        return predicted_map_postprocessed
    
def get_strategy(input_dictionary):
        # prepare GPU id list
    gpu_ids = input_dictionary["gpu_ids"]
    verbose = input_dictionary["verbose"]
    if gpu_ids is None:
        print("No GPU id specified, running on CPU")
        print("If you want to use GPUs, please specify the GPU id(s) using the --gpu_ids flag")
        print("This may take a while...")
        cuda_visible_devices_string = ""
        mirrored_strategy = "cpu"
    else:
        cuda_visible_devices_string = ",".join([str(gpu_id) for gpu_id in gpu_ids])
        #print("Setting CUDA_VISIBLE_DEVICES to {}".format(os.environ["CUDA_VISIBLE_DEVICES"]))
        if verbose:
            print("\tGPU ids: {}".format(cuda_visible_devices_string))
        mirrored_strategy = tf.distribute.MirroredStrategy()
    
    return mirrored_strategy, cuda_visible_devices_string
    
    
def load_emmernet_model(input_dictionary):
    import os
    ## Ignore DeprecationWarning
    import warnings
    warnings.filterwarnings("ignore")
    emmernet_type = input_dictionary["trained_model"]
    emmernet_model_folder = input_dictionary["emmernet_model_folder"]
    model_path = input_dictionary["model_path"]
    verbose = input_dictionary["verbose"]
    
    os.environ["CUDA_VISIBLE_DEVICES"] = input_dictionary["cuda_visible_devices_string"]
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
        emmernet_model_path = os.path.join(emmernet_folder_path, "emmernet", "EMmerNet_MBfa.hdf5")
    elif emmernet_type == "model_free":
        emmernet_model_path = os.path.join(emmernet_folder_path,"emmernet", "EMmerNet_MFfa.hdf5")
    elif emmernet_type == "ensemble":
        emmernet_model_path = os.path.join(emmernet_folder_path, "emmernet","EMmerNet_MBMF.hdf5")
    elif emmernet_type == "hybrid":
        emmernet_model_path = os.path.join(emmernet_folder_path,"emmernet", "epsilon_hybrid_model_4_final_epoch_15.hdf5")
    elif emmernet_type == "model_based_no_freqaug":
        emmernet_model_path = os.path.join(emmernet_folder_path,"emmernet", "EMmerNet_MB.hdf5")
    elif emmernet_type == "emmernet_high_context":
        emmernet_model_path = os.path.join(emmernet_folder_path,"emmernet", "EMmerNet_highContext.hdf5")
    elif emmernet_type == "emmernet_low_context":
        emmernet_model_path = os.path.join(emmernet_folder_path,"emmernet", "EMmerNet_lowContext.hdf5")
    else:
        raise ValueError("Invalid emmernet_type")
    
    if model_path is not None:
        emmernet_model = load_model(model_path, custom_objects={
                                'GroupNormalization': GroupNormalization, \
                                'reducePhysicsBasedLoss': reducePhysicsBasedLoss,
                                'PhysicsBasedMetric': PhysicsBasedMetric,
                                'DataBasedMetric': DataBasedMetric})
        
    else:
        emmernet_model = load_model(emmernet_model_path)
        
    if verbose:
        if model_path is None:
            print("\tEMmerNet model loaded: {}".format(emmernet_type))
        else:
            print("\tEMmerNet model loaded from: {}".format(model_path))
    
    return emmernet_model

    
def run_emmernet_batch(input_dictionary, emmernet_model, mirrored_strategy):
    # collect inputs from input dictionary
    monte_carlo = input_dictionary["monte_carlo"]
    monte_carlo_iterations = input_dictionary["monte_carlo_iterations"]
    batch_size = input_dictionary["batch_size"]
    cubes = input_dictionary["cubes_array"]
    cuda_visible_devices_string = input_dictionary["cuda_visible_devices_string"]
    physics_based = input_dictionary["physics_based"]
    print("Running EMmerNet on {} cubes".format(len(cubes)))
    input_dictionary["logger"].info("Running EMmerNet on {} cubes".format(len(cubes)))
    input_dictionary["logger"].info("Mirrored strategy: {}".format(mirrored_strategy))
    input_dictionary["logger"].info("CUDA_VISIBLE_DEVICES: {}".format(cuda_visible_devices_string))
    
    with RedirectStdoutToLogger(input_dictionary["logger"], show_progress=False):
        emmernet_model.summary()
    
    if mirrored_strategy == "cpu":
        if monte_carlo:
            cubes_predicted_mean, cubes_predicted_var, cubes_predicted_total = run_emmernet_cpu_monte_carlo(cubes, emmernet_model, batch_size, monte_carlo_iterations)
            cubes_predicted_mean = np.squeeze(cubes_predicted_mean, axis=-1)
            cubes_predicted_var = np.squeeze(cubes_predicted_var, axis=-1)
            cubes_predicted_total = np.squeeze(cubes_predicted_total, axis=-1)
        else:
            cubes_predicted_mean = run_emmernet_cpu(cubes, emmernet_model, batch_size)
            cubes_predicted_var = None
            cubes_predicted_total = None
    elif physics_based:        
        cubes_predicted_potential, cubes_predicted_cd = run_emmernet_batch_physics_based(cubes, emmernet_model, batch_size, mirrored_strategy, cuda_visible_devices_string)
        cubes_predicted_mean = cubes_predicted_potential
        cubes_predicted_var = cubes_predicted_cd
        cubes_predicted_total = cubes_predicted_potential
    else:     
        if monte_carlo:
            
            cubes_predicted_mean, cubes_predicted_var, cubes_predicted_total = run_emmernet_batch_monte_carlo(cubes, emmernet_model, batch_size, monte_carlo_iterations, mirrored_strategy, cuda_visible_devices_string)
            cubes_predicted_mean = np.squeeze(cubes_predicted_mean, axis=-1)
            cubes_predicted_var = np.squeeze(cubes_predicted_var, axis=-1)
            cubes_predicted_total = np.squeeze(cubes_predicted_total, axis=-1)
        else:
            cubes_predicted_mean = run_emmernet_batch_no_monte_carlo(cubes, emmernet_model, batch_size, mirrored_strategy, cuda_visible_devices_string)
            cubes_predicted_mean = np.squeeze(cubes_predicted_mean, axis=-1)
            cubes_predicted_var = None
            cubes_predicted_total = None
            
    input_dictionary["cubes_predicted_mean"] = cubes_predicted_mean
    input_dictionary["cubes_predicted_var"] = cubes_predicted_var
    input_dictionary["cubes_predicted_total"] = cubes_predicted_total
    
    return input_dictionary

def run_emmernet_batch_physics_based(cubes, emmernet_model, batch_size, mirrored_strategy, cuda_visible_devices_string):
    import os 
    import warnings
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", category=DeprecationWarning)
        import tensorflow_datasets as tfds
        import atexit
    
    from tqdm import tqdm
    import sys
    
    tfds.disable_progress_bar()
    
    os.environ["CUDA_VISIBLE_DEVICES"] = cuda_visible_devices_string
    
    cube_size = cubes[0].shape[0]
    cubes_predicted_potential = np.empty((0, cube_size, cube_size, cube_size, 1))
    cubes_predicted_cd = np.empty((0, cube_size, cube_size, cube_size, 1))
    cubes_x = np.expand_dims(cubes, axis=4)
    
    with mirrored_strategy.scope():
        for i in tqdm(np.arange(0,len(cubes),batch_size),desc="Running EMmerNet PB",file=sys.stdout):
            if i+batch_size > len(cubes):
                #i = len(cubes)-batch_size-1 # make sure the last batch is of size batch_size
                batch_size = len(cubes)-i
                assert batch_size > 0, "Batch size is less than 0"
                assert batch_size + i == len(cubes), "Batch size and i do not add up to the number of cubes"
            
            cubes_batch_X = np.empty((batch_size, cube_size, cube_size, cube_size, 1))
            cubes_batch_X = cubes_x[i:i+batch_size,:,:,:,:]
            cubes_batch_predicted_potential_cd = emmernet_model(cubes_batch_X, training=False)
            # split potential and cd based on the last dimension
            cubes_batch_predicted_potential = tf.split(cubes_batch_predicted_potential_cd, 2, axis=-1)[0]
            cubes_batch_predicted_cd = tf.split(cubes_batch_predicted_potential_cd, 2, axis=-1)[1]
            cubes_batch_predicted_potential_numpy = cubes_batch_predicted_potential.numpy()
            cubes_batch_predicted_cd_numpy = cubes_batch_predicted_cd.numpy()
            
            cubes_predicted_potential = np.append(cubes_predicted_potential, cubes_batch_predicted_potential_numpy, axis=0)
            cubes_predicted_cd = np.append(cubes_predicted_cd, cubes_batch_predicted_cd_numpy, axis=0)
    
    atexit.register(mirrored_strategy._extended._collective_ops._pool.close)
    
    return cubes_predicted_potential, cubes_predicted_cd
        
def run_emmernet_batch_monte_carlo(cubes, emmernet_model, batch_size, monte_carlo_iterations, mirrored_strategy, cuda_visible_devices_string):
    import os 
    import sys
    import warnings
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", category=DeprecationWarning)
        import tensorflow_datasets as tfds
        import atexit
    
    from tqdm import tqdm
    
    tfds.disable_progress_bar()
    
    os.environ["CUDA_VISIBLE_DEVICES"] = cuda_visible_devices_string
    
    cube_size = cubes[0].shape[0]
    cubes_predicted_full_network = np.empty((0, cube_size, cube_size, cube_size, 1))
    cubes_predicted_mean = np.empty((0, cube_size, cube_size, cube_size, 1))
    cubes_predicted_var = np.empty((0, cube_size, cube_size, cube_size, 1))
    cubes_x = np.expand_dims(cubes, axis=4)
    with mirrored_strategy.scope():
        for i in tqdm(np.arange(0,len(cubes),batch_size),desc="Running MC-EMmerNet", ascii=True, file=sys.stdout):
            if i+batch_size > len(cubes):
                #i = len(cubes)-batch_size-1 # make sure the last batch is of size batch_size
                batch_size = len(cubes)-i
                
                assert batch_size > 0, "Batch size is less than 0"
                assert batch_size + i == len(cubes), "Batch size and i do not add up to the number of cubes"
            
            cubes_batch_X = np.empty((batch_size, cube_size, cube_size, cube_size, 1))
            cubes_batch_X = cubes_x[i:i+batch_size,:,:,:,:]
            cubes_batch_predicted_list = [emmernet_model(cubes_batch_X, training=True) for _ in range(monte_carlo_iterations)]
            # predict again without dropout 
            cubes_batch_predicted = emmernet_model(cubes_batch_X, training=False)
            
            # cubes_batch_predicted_list = [split_potential(cube) for cube in cubes_batch_predicted_list]
            cubes_batch_predicted_numpy = [cube.numpy() for cube in cubes_batch_predicted_list]
            cubes_batch_predicted_mean = np.mean(cubes_batch_predicted_numpy, axis=0)
            cubes_batch_predicted_var = np.var(cubes_batch_predicted_numpy, axis=0)
            cubes_batch_full_network_numpy = cubes_batch_predicted.numpy()

            cubes_predicted_mean = np.append(cubes_predicted_mean, cubes_batch_predicted_mean, axis=0)
            cubes_predicted_var = np.append(cubes_predicted_var, cubes_batch_predicted_var, axis=0)
            cubes_predicted_full_network = np.append(cubes_predicted_full_network, cubes_batch_full_network_numpy, axis=0)
    
    atexit.register(mirrored_strategy._extended._collective_ops._pool.close)
    
    return cubes_predicted_mean, cubes_predicted_var, cubes_predicted_full_network

def compute_mle_mean_variance(data):
    # Reshape data to treat each voxel's values across MC samples as distinct rows
    reshaped_data = data.reshape(data.shape[0], -1).T
    
    # Compute the MLE mean and standard deviation values
    mle_params = np.array([norm.fit(row) for row in reshaped_data])
    
    # Separate and reshape the results back to the original cube shape
    mean_cube = mle_params[:, 0].reshape(data.shape[1:])
    variance_cube = (mle_params[:, 1] ** 2).reshape(data.shape[1:])
    
    return mean_cube, variance_cube
def run_emmernet_batch_no_monte_carlo(cubes, emmernet_model, batch_size, mirrored_strategy, cuda_visible_devices_string):
    import os
    import warnings
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", category=DeprecationWarning)
        import tensorflow_datasets as tfds
        import atexit
    
    from tqdm import tqdm
    import sys

    tfds.disable_progress_bar()

    os.environ["CUDA_VISIBLE_DEVICES"] = cuda_visible_devices_string
    
    cube_size = cubes[0].shape[0]
    cubes_predicted = np.empty((0, cube_size, cube_size, cube_size, 1))
    cubes_x = np.expand_dims(cubes, axis=4)
    
    with mirrored_strategy.scope():
        for i in tqdm(np.arange(0,len(cubes),batch_size),desc="Running EMmerNet",file=sys.stdout):
            if i+batch_size > len(cubes):
                #i = len(cubes)-batch_size-1 # make sure the last batch is of size batch_size
                batch_size = len(cubes)-i
                
                assert batch_size > 0, "Batch size is less than 0"
                assert batch_size + i == len(cubes), "Batch size and i do not add up to the number of cubes"
            
            cubes_batch_X = np.empty((batch_size, cube_size, cube_size, cube_size, 1))
            cubes_batch_X = cubes_x[i:i+batch_size,:,:,:,:]
            cubes_batch_predicted = emmernet_model(cubes_batch_X, training=True)
            cubes_predicted = np.append(cubes_predicted, cubes_batch_predicted, axis=0)
    
    atexit.register(mirrored_strategy._extended._collective_ops._pool.close)
    
    return cubes_predicted
    
    
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

class reducePhysicsBasedLoss(tf.keras.losses.Loss):
        """ custom loss function that reduces physics based loss
        """
        def __init__(self,reduction=tf.keras.losses.Reduction.AUTO, name="reducePhysicsBasedLoss"):
            super().__init__(name="reducePhysicsBasedLoss")
        
        def laplacian_tf(self, tensor, dx=1., dy=1., dz=1.):
            tensor_shape = tf.shape(tensor)
            tensor = tf.reshape(tensor, [-1, tensor_shape[1], tensor_shape[2], tensor_shape[3]])
            
            z, y, x = tf.meshgrid(tf.range(32), tf.range(32), tf.range(32), indexing='ij')
            x = tf.cast(x, dtype=tf.float32) * dx
            y = tf.cast(y, dtype=tf.float32) * dy
            z = tf.cast(z, dtype=tf.float32) * dz

            with tf.GradientTape() as tape2:
                tape2.watch([x, y, z])
                with tf.GradientTape() as tape1:
                    tape1.watch([x, y, z])
                    grad_x = tape1.gradient(tensor, x, unconnected_gradients='zero')
                    grad_y = tape1.gradient(tensor, y, unconnected_gradients='zero')
                    grad_z = tape1.gradient(tensor, z, unconnected_gradients='zero')
                laplacian_x = tape2.gradient(grad_x, x, unconnected_gradients='zero')
                laplacian_y = tape2.gradient(grad_y, y, unconnected_gradients='zero')
                laplacian_z = tape2.gradient(grad_z, z, unconnected_gradients='zero')
                
            laplacian = laplacian_x + laplacian_y + laplacian_z
            laplacian = tf.reshape(laplacian, tensor_shape)
            return laplacian
                
        # def laplacian_tf(self, tensor):
        #     laplace_kernel = tf.constant([[[0, 0, 0], [0, 1, 0], [0, 0, 0]],
        #                                 [[0, 1, 0], [1, -6, 1], [0, 1, 0]],
        #                                 [[0, 0, 0], [0, 1, 0], [0, 0, 0]]], dtype=tf.float32)

        #     laplace_kernel = tf.reshape(laplace_kernel, [3, 3, 3, 1, 1])
        #     tensor = tf.expand_dims(tensor, -1)
        #     laplacian = tf.nn.conv3d(tensor, laplace_kernel, [1, 1, 1, 1, 1], "SAME")
        #     laplacian = tf.squeeze(laplacian, -1)

        #     return laplacian

        # Then, in your loss function, you can use this layer to compute the Laplacian:
        def physics_based_loss(self, y_pred, y_true):
            potential_tf, charge_density_tf = tf.split(y_pred, num_or_size_splits=2, axis=-1)
            #laplacian_layer = LaplacianLayer()
            laplacian_potential_tf = -1 * self.laplacian_tf(potential_tf)
        
            return tf.reduce_mean(tf.square(potential_tf - charge_density_tf))

        def simplified_loss(self, y_true, y_pred):
            print("SHAPE OF Y_PRED: ", y_pred.shape)
            print("SHAPE OF Y_TRUE: ", y_true.shape)
            return tf.reduce_mean(tf.square(y_true - y_pred))
                
        def __call__(self, y_true, y_pred, sample_weight=None):
            return self.physics_based_loss(y_pred=y_pred, y_true=y_true)

class PhysicsBasedMetric(tf.keras.metrics.Metric):
    def __init__(self, name='PhysicsBasedLoss', **kwargs):
        super(PhysicsBasedMetric, self).__init__(name=name, **kwargs)
        self.physics_based_loss = self.add_weight(name='pb_loss', initializer='zeros')

    def update_state(self, y_true, y_pred, sample_weight=None):
        physics_loss = reducePhysicsBasedLoss().physics_based_loss(y_pred, y_true)
        self.physics_based_loss.assign(physics_loss)

    def result(self):
        return self.physics_based_loss

class DataBasedMetric(tf.keras.metrics.Metric):
    def __init__(self, name='DataBasedLoss', **kwargs):
        super(DataBasedMetric, self).__init__(name=name, **kwargs)
        self.data_based_loss = self.add_weight(name='db_loss', initializer='zeros')

    def update_state(self, y_true, y_pred, sample_weight=None):
        data_loss = reducePhysicsBasedLoss().data_based_loss(y_pred, y_true)
        self.data_based_loss.assign(data_loss)

    def result(self):
        return self.data_based_loss



def run_emmernet_cpu(cubes, emmernet_model, batch_size):
    ## Run the model on the cube
    ## Ignore DeprecationWarning
    import warnings
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", category=DeprecationWarning)
        import tensorflow_datasets as tfds
        import atexit
    from tqdm import tqdm
    import os
    import sys

    print("Running EMmerNet on {} cubes using CPU".format(len(cubes)))

    tfds.disable_progress_bar()
    
    cube_size = cubes[0].shape[0]
    cubes = np.array(cubes)
    cubes_x = np.expand_dims(cubes, axis=4)
    
    cubes_predicted = np.empty((0, cube_size, cube_size, cube_size, 1))

    for i in tqdm(np.arange(0,len(cubes),batch_size),desc="Running EMmerNet", file=sys.stdout):
        if i+batch_size > len(cubes):
            #i = len(cubes)-batch_size-1 # make sure the last batch is of size batch_size
            batch_size = len(cubes)-i
            
            assert batch_size > 0, "Batch size is less than 0"
            assert batch_size + i == len(cubes), "Batch size and i do not add up to the number of cubes"
        
        cubes_batch_X = np.empty((batch_size, cube_size, cube_size, cube_size, 1))
        cubes_batch_X = cubes_x[i:i+batch_size,:,:,:,:]

        cubes_batch_predicted = emmernet_model.predict(x=cubes_batch_X, batch_size=batch_size, verbose=0)
        cubes_predicted = np.append(cubes_predicted, cubes_batch_predicted, axis=0)
    
    cubes_predicted = np.squeeze(cubes_predicted, axis=-1)
    return cubes_predicted

def run_emmernet_cpu_monte_carlo(cubes, emmernet_model, batch_size, monte_carlo_iterations):
    import os 
    import warnings
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", category=DeprecationWarning)
        import tensorflow_datasets as tfds
        import atexit
    
    from tqdm import tqdm
    import sys

    tfds.disable_progress_bar()
    
    
    cube_size = cubes[0].shape[0]
    cubes_predicted_full_network = np.empty((0, cube_size, cube_size, cube_size, 1))
    cubes_predicted_mean = np.empty((0, cube_size, cube_size, cube_size, 1))
    cubes_predicted_var = np.empty((0, cube_size, cube_size, cube_size, 1))
    cubes_x = np.expand_dims(cubes, axis=4)

    for i in tqdm(np.arange(0,len(cubes),batch_size),desc="Running EMmerNet", file=sys.stdout):
        if i+batch_size > len(cubes):
            #i = len(cubes)-batch_size-1 # make sure the last batch is of size batch_size
            batch_size = len(cubes)-i
            
            assert batch_size > 0, "Batch size is less than 0"
            assert batch_size + i == len(cubes), "Batch size and i do not add up to the number of cubes"
        
        cubes_batch_X = np.empty((batch_size, cube_size, cube_size, cube_size, 1))
        cubes_batch_X = cubes_x[i:i+batch_size,:,:,:,:]
        cubes_batch_predicted_list = [emmernet_model(cubes_batch_X, training=True) for _ in range(monte_carlo_iterations)]
        # predict again without dropout 
        cubes_batch_predicted = emmernet_model(cubes_batch_X, training=False)
        
        # cubes_batch_predicted_list = [split_potential(cube) for cube in cubes_batch_predicted_list]
        cubes_batch_predicted_numpy = [cube.numpy() for cube in cubes_batch_predicted_list]
        cubes_batch_predicted_mean = np.mean(cubes_batch_predicted_numpy, axis=0)
        cubes_batch_predicted_var = np.var(cubes_batch_predicted_numpy, axis=0)
        cubes_batch_full_network_numpy = cubes_batch_predicted.numpy()

        cubes_predicted_mean = np.append(cubes_predicted_mean, cubes_batch_predicted_mean, axis=0)
        cubes_predicted_var = np.append(cubes_predicted_var, cubes_batch_predicted_var, axis=0)
        cubes_predicted_full_network = np.append(cubes_predicted_full_network, cubes_batch_full_network_numpy, axis=0)

    return cubes_predicted_mean, cubes_predicted_var, cubes_predicted_full_network
    
