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
from locscale.emmernet.emmernet_functions_new import standardize_map, minmax_normalize_map, get_cubes, assemble_cubes, replace_cubes_in_dictionary,\
                                                    load_smoothened_mask, show_signal_cubes
                                                    
from locscale.emmernet.utils import compute_local_phase_correlations, plot_phase_correlations
import tensorflow as tf
import numpy as np
import os
from tqdm import tqdm

os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3'  
tf.random.set_seed(42)

def run_emmernet(input_dictionary):
    ## Ignore DeprecationWarning
    import warnings
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", category=DeprecationWarning)
        import tensorflow as tf
        from tensorflow.keras.models import load_model
    
    EMMERNET_CUBE_SIZE=input_dictionary["cube_size"]
    
    ## Get the input map path
    emmap_path = input_dictionary["emmap_path"]
    processing_files_folder = os.path.dirname(emmap_path)
    mask_path = input_dictionary["xyz_mask_path"]
    emmernet_type = input_dictionary["trained_model"]
    stride = input_dictionary["stride"]
    batch_size = input_dictionary["batch_size"]
    gpu_ids = input_dictionary["gpu_ids"]
    verbose = input_dictionary["verbose"]
    emmernet_model_folder = input_dictionary["emmernet_model_folder"]
    target_map_path = input_dictionary["target_map_path"]
    model_path = input_dictionary["model_path"]
    monte_carlo = input_dictionary["monte_carlo"]
    monte_carlo_iterations = input_dictionary["monte_carlo_iterations"]

    if target_map_path is not None:
        target_map, _ = load_map(target_map_path)
        target_map_present = True
    else:
        target_map_present = False
    

    emmap, apix = load_map(emmap_path)
    mask, _ = load_smoothened_mask(mask_path)
    
    input_map_shape = emmap.shape
    if verbose:
        print("Emmap loaded from: {}".format(emmap_path))
        print("Emmap shape: {}".format(emmap.shape))
        print("Pixelsize read as: {:.2f}".format(apix))

        print("1) Pre-processing commencing...")
    
    ## Preprocess

    emmap_preprocessed = preprocess_map(emmap, apix)
    mask_preprocessed = preprocess_map(mask, apix, standardize=False)
    if target_map_present:
        target_map_preprocessed = preprocess_map(target_map, apix)
    if verbose:
        print("\tPreprocessing complete")
        print("\tPre-processed map shape: {}".format(emmap_preprocessed.shape))
        print("2) Prediction commencing...")

    cubes_dictionary, cubes_array, signal_cubes = get_cubes(emmap_preprocessed, cube_size=EMMERNET_CUBE_SIZE, step_size=stride, mask=mask_preprocessed)
    show_signal_cubes(signal_cubes, emmap_preprocessed.shape, save_path=os.path.join(processing_files_folder, "signal_cubes_resampled.mrc"), apix=1)

    if target_map_present:
        target_cubes_dictionary, target_cubes_array = get_cubes(target_map_preprocessed, cube_size=EMMERNET_CUBE_SIZE, step_size=stride, mask=mask_preprocessed)
        
    if verbose:
        print("\tCubes extracted")
        print("\tNumber of cubes: {}".format(len(cubes_array)))
    ## Load the model
    
    # prepare GPU id list
    if gpu_ids is None:
        print("No GPU id specified, running on CPU")
        print("If you want to use GPUs, please specify the GPU id(s) using the --gpu_ids flag")
        print("This may take a while...")
        # set device to CPU
        os.environ["CUDA_VISIBLE_DEVICES"] = ""
        device = "/cpu:0"
        print("Setting CUDA_VISIBLE_DEVICES to {}".format(os.environ["CUDA_VISIBLE_DEVICES"]))
        print("Using device: {}".format(device))
        mirrored_strategy = "cpu"
    else:
        os.environ["CUDA_VISIBLE_DEVICES"] = ",".join([str(gpu_id) for gpu_id in gpu_ids])
        print("Setting CUDA_VISIBLE_DEVICES to {}".format(os.environ["CUDA_VISIBLE_DEVICES"]))
        gpu_id_list = ["/gpu:"+str(gpu_id) for gpu_id in gpu_ids]
        if verbose:
            print("\tGPU ids: {}".format(gpu_id_list))
        
        mirrored_strategy = tf.distribute.MirroredStrategy()

    emmernet_model = load_emmernet_model(emmernet_type, emmernet_model_folder, model_path)
    if verbose:
        if model_path is None:
            print("\tEMmerNet model loaded: {}".format(emmernet_type))
        else:
            print("\tEMmerNet model loaded from: {}".format(model_path))

    ## Run EMmerNet using GPUs
    

    
    if monte_carlo:
        predicted_cubes_mean, predicted_cubes_var = run_emmernet_batch(cubes_array, emmernet_model, mirrored_strategy, \
                                                        batch_size=batch_size, monte_carlo_iterations=monte_carlo_iterations)
        predicted_cubes_mean_dictionary = replace_cubes_in_dictionary(predicted_cubes_mean, cubes_dictionary)
        predicted_cubes_var_dictionary = replace_cubes_in_dictionary(predicted_cubes_var, cubes_dictionary)
        length_of_predicted_cubes = len(predicted_cubes_mean)
    else:
        predicted_cubes_potential = run_emmernet_batch(cubes_array, emmernet_model, mirrored_strategy, batch_size=batch_size)
        predicted_cubes_dictionary = replace_cubes_in_dictionary(predicted_cubes_potential, cubes_dictionary)
        # predicted_cubes_charge_density_dictionary = replace_cubes_in_dictionary(predicted_cubes_charge_density, cubes_dictionary)
        length_of_predicted_cubes = len(predicted_cubes_potential)
    if target_map_present:
        # Compute phase correlations between predicted and target cubes
        
        phase_correlations, freq = compute_local_phase_correlations(target_cubes=target_cubes, predicted_cubes=predicted_cubes, apix=apix, temp_folder=processing_files_folder)
        phase_correlations_fig_path = os.path.join(processing_files_folder, "phase_correlations.png")
        fig = plot_phase_correlations(phase_correlations, freq)
        fig.savefig(phase_correlations_fig_path)
        
    if verbose:
        print("\tEMmerNet prediction complete")
        print("\tNumber of predicted cubes: {}".format(length_of_predicted_cubes))
    ## Merge the predicted cubes sequence
    
    if monte_carlo:
        predicted_map_mean = assemble_cubes(predicted_cubes_mean_dictionary, emmap_preprocessed.shape[0])
        predicted_map_var = assemble_cubes(predicted_cubes_var_dictionary, emmap_preprocessed.shape[0])
        final_shape_prediction = predicted_map_mean.shape
    else:
        predicted_map_potential = assemble_cubes(predicted_cubes_dictionary, emmap_preprocessed.shape[0])
        #predicted_map_charge_density = assemble_cubes(predicted_cubes_charge_density_dictionary, emmap_preprocessed.shape[0])
        final_shape_prediction = predicted_map_potential.shape
    
    if verbose:
        print("\tPredicted map assembled")
        print("\tPredicted map shape: {}".format(final_shape_prediction))
        print("3) Post-processing commencing...")
    
    
    ## Postprocess
    if monte_carlo:
        predicted_map_mean_postprocessed = postprocess_map(predicted_map_mean, apix, output_shape=input_map_shape)
        predicted_map_var_postprocessed = postprocess_map(predicted_map_var, apix, output_shape=input_map_shape)
        post_processed_map_shape = predicted_map_mean_postprocessed.shape
    else:
        predicted_map_postprocessed = postprocess_map(predicted_map_potential, apix, output_shape=input_map_shape)
        #predicted_map_charge_density_postprocessed = postprocess_map(predicted_map_charge_density, apix, output_shape=input_map_shape)
        post_processed_map_shape = predicted_map_postprocessed.shape
        
    if verbose:
        print("\tPost-processing complete")
        print("\tPost-processed map shape: {}".format(post_processed_map_shape))

    #return predicted_map_postprocessed

    if monte_carlo:
        emmernet_output_dictionary = {"output_mean":predicted_map_mean_postprocessed, "output_var":predicted_map_var_postprocessed}
    else:
        emmernet_output_dictionary = {"output":predicted_map_postprocessed}#, "output_charge_density":predicted_map_charge_density_postprocessed}

    return emmernet_output_dictionary


def load_emmernet_model(emmernet_type, emmernet_model_folder=None, model_path=None):
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
        emmernet_model_path = os.path.join(emmernet_folder_path, "EMmerNet_MBfa.hdf5")
    elif emmernet_type == "model_free":
        emmernet_model_path = os.path.join(emmernet_folder_path, "EMmerNet_MFfa.hdf5")
    elif emmernet_type == "ensemble":
        emmernet_model_path = os.path.join(emmernet_folder_path, "EMmerNet_MBMF.hdf5")
    elif emmernet_type == "hybrid":
        emmernet_model_path = os.path.join(emmernet_folder_path, "epsilon_hybrid_model_4_final_epoch_15.hdf5")
    elif emmernet_type == "model_based_no_freqaug":
        emmernet_model_path = os.path.join(emmernet_folder_path, "EMmerNet_MB.hdf5")
    elif emmernet_type == "hybrid_model_map_target":
        emmernet_model_path = os.path.join(emmernet_folder_path, "hybrid_model_map_target_noaugmentation_30k_bs2_dropout_epoch_15.hdf5")
    else:
        raise ValueError("Invalid emmernet_type")
    
    if model_path is not None:
        emmernet_model = load_model(model_path, custom_objects={
                                'GroupNormalization': GroupNormalization, \
                                'reducePhysicsBasedLoss': reducePhysicsBasedLoss,
                                'PhysicsBasedMetric': PhysicsBasedMetric,
                                'DataBasedMetric': DataBasedMetric})
        return emmernet_model
    else:
        emmernet_model = load_model(emmernet_model_path)
        return emmernet_model

def run_emmernet_cpu(cubes, emmernet_model, batch_size, monte_carlo_iterations=None):
    ## Run the model on the cube
    ## Ignore DeprecationWarning
    import warnings
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", category=DeprecationWarning)
        import tensorflow_datasets as tfds
        import atexit
    from tqdm import tqdm
    import os
    
    print("Running EMmerNet on {} cubes using CPU".format(len(cubes)))
    monte_carlo = monte_carlo_iterations is not None

    tfds.disable_progress_bar()
    cube_size = cubes[0].shape[0]
    cubes = np.array(cubes)
    cubes_x = np.expand_dims(cubes, axis=4)
    if monte_carlo:
        cubes_predicted_mean = np.empty((0, cube_size, cube_size, cube_size, 1))
        cubes_predicted_var = np.empty((0, cube_size, cube_size, cube_size, 1))
    else:
        cubes_predicted = np.empty((0, cube_size, cube_size, cube_size, 1))
    mean_values = []
    var_values = []
    
    split_potential = lambda x: tf.split(x, num_or_size_splits=2, axis=-1)[0]
    split_charge_density = lambda x: tf.split(x, num_or_size_splits=2, axis=-1)[1]
    cubes_charge_density = np.empty((0, cube_size, cube_size, cube_size, 1))
    
    
    for i in tqdm(np.arange(0,len(cubes),batch_size),desc="Running EMmerNet"):
        if i+batch_size > len(cubes):
            #i = len(cubes)-batch_size-1 # make sure the last batch is of size batch_size
            batch_size = len(cubes)-i
            
            assert batch_size > 0, "Batch size is less than 0"
            assert batch_size + i == len(cubes), "Batch size and i do not add up to the number of cubes"
        
        cubes_batch_X = np.empty((batch_size, cube_size, cube_size, cube_size, 1))
        cubes_batch_X = cubes_x[i:i+batch_size,:,:,:,:]

        if monte_carlo:
            cubes_batch_predicted_list = [emmernet_model(cubes_batch_X, training=True) for _ in range(monte_carlo_iterations)]
            # cubes_batch_predicted_list = [split_potential(cube) for cube in cubes_batch_predicted_list]
            cubes_batch_predicted_numpy = [cube.numpy() for cube in cubes_batch_predicted_list]
            cubes_batch_predicted_mean = np.mean(cubes_batch_predicted_numpy, axis=0)
            cubes_batch_predicted_var = np.var(cubes_batch_predicted_numpy, axis=0)

            mean_values.append(cubes_batch_predicted_mean.mean())
            var_values.append(cubes_batch_predicted_var.mean())
            cubes_predicted_mean = np.append(cubes_predicted_mean, cubes_batch_predicted_mean, axis=0)
            cubes_predicted_var = np.append(cubes_predicted_var, cubes_batch_predicted_var, axis=0)
        else:
            cubes_batch_predicted = emmernet_model.predict(x=cubes_batch_X, batch_size=batch_size, verbose=0)
            
            cubes_batch_predicted_first_channel = split_potential(cubes_batch_predicted)
            cubes_batch_second_channel = split_charge_density(cubes_batch_predicted)
            cubes_predicted = np.append(cubes_predicted, cubes_batch_predicted_first_channel, axis=0)
            cubes_charge_density = np.append(cubes_charge_density, cubes_batch_second_channel, axis=0)
    
    # squeeze cubes to 3 dimensions
    if monte_carlo:
        cubes_predicted_mean = np.squeeze(cubes_predicted_mean, axis=-1)
        cubes_predicted_var = np.squeeze(cubes_predicted_var, axis=-1)
        return cubes_predicted_mean, cubes_predicted_var
    else:
        cubes_predicted = np.squeeze(cubes_predicted, axis=-1)
        cubes_charge_density = np.squeeze(cubes_charge_density, axis=-1)
        return cubes_predicted, cubes_charge_density
    
    
    
def run_emmernet_batch(cubes, emmernet_model, mirrored_strategy, batch_size, monte_carlo_iterations=None):
    ## Run the model on the cube
    ## Ignore DeprecationWarning
    
    if mirrored_strategy == "cpu":
        return run_emmernet_cpu(cubes, emmernet_model, batch_size, monte_carlo_iterations)
    
    import warnings
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", category=DeprecationWarning)
        import tensorflow_datasets as tfds
        import atexit
    from tqdm import tqdm
    import os
    
    print("Running EMmerNet on {} cubes".format(len(cubes)))
    monte_carlo = monte_carlo_iterations is not None

    tfds.disable_progress_bar()
    cube_size = cubes[0].shape[0]
    cubes = np.array(cubes)
    cubes_x = np.expand_dims(cubes, axis=4)
    if monte_carlo:
        cubes_predicted_mean = np.empty((0, cube_size, cube_size, cube_size, 1))
        cubes_predicted_var = np.empty((0, cube_size, cube_size, cube_size, 1))
    else:
        cubes_predicted = np.empty((0, cube_size, cube_size, cube_size, 1))
    mean_values = []
    var_values = []
    
    split_potential = lambda x: tf.split(x, num_or_size_splits=2, axis=-1)[0]
    split_charge_density = lambda x: tf.split(x, num_or_size_splits=2, axis=-1)[1]
    cubes_charge_density = np.empty((0, cube_size, cube_size, cube_size, 1))
    
    
    with mirrored_strategy.scope():
        for i in tqdm(np.arange(0,len(cubes),batch_size),desc="Running EMmerNet"):
            if i+batch_size > len(cubes):
                #i = len(cubes)-batch_size-1 # make sure the last batch is of size batch_size
                batch_size = len(cubes)-i
                
                assert batch_size > 0, "Batch size is less than 0"
                assert batch_size + i == len(cubes), "Batch size and i do not add up to the number of cubes"
            
            cubes_batch_X = np.empty((batch_size, cube_size, cube_size, cube_size, 1))
            cubes_batch_X = cubes_x[i:i+batch_size,:,:,:,:]

            if monte_carlo:
                cubes_batch_predicted_list = [emmernet_model(cubes_batch_X, training=True) for _ in range(monte_carlo_iterations)]
                # cubes_batch_predicted_list = [split_potential(cube) for cube in cubes_batch_predicted_list]
                cubes_batch_predicted_numpy = [cube.numpy() for cube in cubes_batch_predicted_list]
                cubes_batch_predicted_mean = np.mean(cubes_batch_predicted_numpy, axis=0)
                cubes_batch_predicted_var = np.var(cubes_batch_predicted_numpy, axis=0)
    
                mean_values.append(cubes_batch_predicted_mean.mean())
                var_values.append(cubes_batch_predicted_var.mean())
                cubes_predicted_mean = np.append(cubes_predicted_mean, cubes_batch_predicted_mean, axis=0)
                cubes_predicted_var = np.append(cubes_predicted_var, cubes_batch_predicted_var, axis=0)
            else:
                cubes_batch_predicted = emmernet_model.predict(x=cubes_batch_X, batch_size=batch_size, verbose=0)
                
                #cubes_batch_predicted_first_channel = split_potential(cubes_batch_predicted)
                #cubes_batch_second_channel = split_charge_density(cubes_batch_predicted)
                cubes_predicted = np.append(cubes_predicted, cubes_batch_predicted, axis=0)
                #cubes_charge_density = np.append(cubes_charge_density, cubes_batch_second_channel, axis=0)
        
        
    # close the mirrored strategy's multiprocessing ThreadPool explicitly
    atexit.register(mirrored_strategy._extended._collective_ops._pool.close)

    # squeeze cubes to 3 dimensions
    if monte_carlo:
        cubes_predicted_mean = np.squeeze(cubes_predicted_mean, axis=-1)
        cubes_predicted_var = np.squeeze(cubes_predicted_var, axis=-1)
        return cubes_predicted_mean, cubes_predicted_var
    else:
        cubes_predicted = np.squeeze(cubes_predicted, axis=-1)
        #cubes_charge_density = np.squeeze(cubes_charge_density, axis=-1)
        return cubes_predicted#, cubes_charge_density


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
    ## MinMax normalize the map
    #predicted_map_normalized = minmax_normalize_map(predicted_map_resampled)
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
