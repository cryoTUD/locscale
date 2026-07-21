#
# Delft University of Technology (TU Delft) hereby disclaims all copyright interest in the program 'LocScale'
# written by the Author(s).
# Copyright (C) 2021 Alok Bharadwaj and Arjen J. Jakobi
# This software may be modified and distributed under the terms of the BSD license. 
# You should have received a copy of the BSD 3-clause license along with this program (see LICENSE file file for details).
# If not see https://opensource.org/license/bsd-3-clause/.
#

import numpy as np
import os
import sys
from locscale.include.emmer.ndimage.map_utils import save_as_mrc

def run_window_function_including_scaling(parsed_inputs_dict):
    """
    This is a function which performs high level data processing for Locscale

    """
    from locscale.utils.general import get_xyz_locs_and_indices_after_edge_cropping_and_masking
    from locscale.utils.general import save_list_as_map, put_scaled_voxels_back_in_original_volume_including_padding
    from locscale.utils.general import merge_sequence_of_sequences, split_sequence_evenly, write_out_final_volume_window_back_if_required
    from joblib import Parallel, delayed
    from locscale.utils.amplitude_scaling import local_amplitude_scaling
    ###############################################################################
    # Stage 1: Collect inputs from the dictionary
    ###############################################################################

    scaling_dictionary = parsed_inputs_dict
    ###############################################################################
    # Stage 2: Extract masked locations and indices from the mask
    ###############################################################################
    
    corner_voxels, masked_indices, map_shape = get_xyz_locs_and_indices_after_edge_cropping_and_masking(
        scaling_dictionary['mask'], scaling_dictionary['wn'])
    
    reference_map = scaling_dictionary["modmap"]
    target_map = scaling_dictionary["emmap"]
    window_size = scaling_dictionary["wn"]
    chunk = scaling_dictionary["chunk"]
    
    # sharpened_vals = local_amplitude_scaling(reference_map, target_map, corner_voxels, window_size=window_size, chunk=chunk)

    masked_xyz_locs_split = split_sequence_evenly(corner_voxels, scaling_dictionary['number_processes'])

    scaling_dictionary["masked_indices"] = masked_indices
    scaling_dictionary["map_shape"] = map_shape
    

    scaling_dictionary_split = {}
    for i in range(scaling_dictionary['number_processes']):
        scaling_dictionary_split[i] = scaling_dictionary.copy()
        scaling_dictionary_split[i]["corner_voxels"] = masked_xyz_locs_split[i]
        scaling_dictionary_split[i]["use_mpi"] = False
    ###############################################################################
    # Stage 3: Run the window function to get sharpened values and bfactor information
    ###############################################################################
    # Use joblib to parallelize the window function 
    if scaling_dictionary['number_processes'] > 1:
        results = Parallel(n_jobs=scaling_dictionary['number_processes'])(
            delayed(local_amplitude_scaling)(
                reference_map, 
                target_map, 
                scaling_dictionary_split[i]["corner_voxels"],
                window_size=window_size,
                chunk=chunk
            ) for i in range(scaling_dictionary['number_processes'])
        )
    else:
        scaling_dictionary_split[0]["use_mpi"] = False
        results = [local_amplitude_scaling(
            reference_map, 
            target_map, 
            scaling_dictionary_split[i]["corner_voxels"],
            window_size=window_size,
            chunk=chunk
            )]
    
    # ###############################################################################
    # # Stage 4: Merge the results from the parallelized window function
    # ###############################################################################
    if scaling_dictionary['number_processes'] > 1:
        sharpened_vals = merge_sequence_of_sequences([results[i] for i in range(scaling_dictionary['number_processes'])])
    else:
        sharpened_vals = results[0]
        

    # ###############################################################################
    # # Stage 5: Put the sharpened values back in the original volume
    # ###############################################################################

    map_scaled = put_scaled_voxels_back_in_original_volume_including_padding(sharpened_vals, masked_indices, map_shape)
    
    ###############################################################################
    # Stage 6: Return the scaled map
    ###############################################################################
    return map_scaled

def run_window_function_including_scaling_mpi(parsed_inputs_dict):
    """
    This is a function which performs high level data processing for Locscale in a MPI environment

    """

    from mpi4py import MPI
    from locscale.utils.general import get_xyz_locs_and_indices_after_edge_cropping_and_masking
    from locscale.utils.general import save_list_as_map, merge_sequence_of_sequences, split_sequence_evenly
    from locscale.utils.general import put_scaled_voxels_back_in_original_volume_including_padding
    from locscale.utils.amplitude_scaling import local_amplitude_scaling
                                   
    ###############################################################################
    # Stage 1: Collect inputs from the dictionary
    ###############################################################################
    scaling_dictionary_mpi = parsed_inputs_dict
    
    ###############################################################################
    # Stage 1a: Setup MPI environment
    ###############################################################################
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    size = comm.Get_size()

    ###############################################################################
    # Stage 2: Extract masked locations and indices from the mask from root node
    ###############################################################################
    if rank == 0:
        masked_xyz_locs, masked_indices, map_shape = \
        get_xyz_locs_and_indices_after_edge_cropping_and_masking(scaling_dictionary_mpi['mask'], scaling_dictionary_mpi['wn'])

        zs, ys, xs = masked_xyz_locs.T
        zs = split_sequence_evenly(zs, size)
        ys = split_sequence_evenly(ys, size)
        xs = split_sequence_evenly(xs, size)
    else:
        zs = None
        ys = None
        xs = None

    zs = comm.scatter(zs, root=0)
    ys = comm.scatter(ys, root=0)
    xs = comm.scatter(xs, root=0)

    masked_xyz_locs = np.column_stack((zs, ys, xs))

    process_name = 'LocScale process {0} of {1}'.format(rank + 1, size)
    scaling_dictionary_mpi['masked_xyz_locs'] = masked_xyz_locs
    scaling_dictionary_mpi["use_mpi"] = True
    if rank == 0:
        scaling_dictionary_mpi['masked_indices'] = masked_indices
        scaling_dictionary_mpi['map_shape'] = map_shape
        
    ###############################################################################
    # Stage 3: Run the window function to get sharpened values and bfactor information
    ###############################################################################

    sharpened_vals = local_amplitude_scaling(
        scaling_dictionary_mpi['modmap'],
        scaling_dictionary_mpi['emmap'],
        scaling_dictionary_mpi['masked_xyz_locs'],
        window_size=scaling_dictionary_mpi['wn'],
        chunk=scaling_dictionary_mpi['chunk']
    )

    
    ###############################################################################
    # Stage 4: Put the sharpened values back into the original volume
    ###############################################################################
    
    ###############################################################################
    # Stage 4a: Gather the computed values from all nodes to the root node
    ###############################################################################
    sharpened_vals = comm.gather(sharpened_vals, root=0)

    if rank == 0:
        sharpened_vals = merge_sequence_of_sequences(sharpened_vals)        
        map_scaled = put_scaled_voxels_back_in_original_volume_including_padding(np.array(sharpened_vals),
        masked_indices, map_shape)

        ###########################################################################
        # Stage 5: Save the processing files 
        ###########################################################################
    else:
        map_scaled = None

    ######## Wait for all processes to finish #########
    comm.barrier()

    ###############################################################################
    # Stage 6: Return the scaled map
    ###############################################################################
    return map_scaled, rank




