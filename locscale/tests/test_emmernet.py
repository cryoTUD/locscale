#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 13 00:01:47 2021

@author: alok
"""

import unittest
import numpy as np
import os

class test_emmernet(unittest.TestCase):
    def setUp(self):
        from locscale.emmernet.utils import check_and_download_emmernet_model

        emmernet_model_folder = check_and_download_emmernet_model(verbose=True)
        self.assertTrue(emmernet_model_folder is not None)
        self.emmernet_model_folder = emmernet_model_folder
        self.inputs_dictionary = {
            "trained_model" : "emmernet_high_context",
            "emmernet_model_folder" : self.emmernet_model_folder,
            "model_path" : None,
            "verbose" : True,
            "cuda_visible_devices_string" : "",
        }
        self.emmap_shape = (252,252,252)
        self.cube_size = 32
        self.stride = 16
        self.emmap = np.random.randint(0,255,self.emmap_shape)
        # get a random sphere in the map
        self.mask = np.zeros_like(self.emmap)
        self.mask[self.emmap_shape[0]//2-25:self.emmap_shape[0]//2+25,
                self.emmap_shape[1]//2-25:self.emmap_shape[1]//2+25,
                self.emmap_shape[2]//2-25:self.emmap_shape[2]//2+25] = 1
        self.masked_emmap = self.emmap * self.mask

        
    def test_map_chunking(self):
        from locscale.emmernet.emmernet_functions import get_cubes, assemble_cubes
        import numpy as np
        print("Testing: map_chunking")

        
        # test 1: get_cubes
        cubes_dictionary, cubes_array, filtered_signal_cubecenters = get_cubes(self.masked_emmap, self.stride, self.cube_size, self.mask)

        # test 2: assemble_cubes
        
        assembled_cube = assemble_cubes(cubes_dictionary, self.emmap_shape[0], average=True)

        # test 3: check if assembled cube is same as original emmap
        self.assertTrue(np.allclose(assembled_cube, self.masked_emmap))
        self.assertTrue(np.squeeze(cubes_array[0]).shape == (self.cube_size, self.cube_size, self.cube_size))
    
    def test_load_emmernet_model(self):
        import warnings
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            from locscale.emmernet.run_emmernet import load_emmernet_model

        inputs_dictionary = {
            "trained_model" : "emmernet_high_context",
            "emmernet_model_folder" : self.emmernet_model_folder,
            "model_path" : None,
            "verbose" : True,
            "cuda_visible_devices_string" : "",
        }
        
        inputs_dictionary_2 = {
            "trained_model" : "emmernet_low_context",
            "emmernet_model_folder" : self.emmernet_model_folder,
            "model_path" : None,
            "verbose" : True,
            "cuda_visible_devices_string" : "",
        }
        
        emmernet_model_1 = load_emmernet_model(inputs_dictionary)
        emmernet_model_2 = load_emmernet_model(inputs_dictionary_2)

        self.assertTrue(emmernet_model_1 is not None)
        self.assertTrue(emmernet_model_2 is not None)
        self.test_input_dict = inputs_dictionary
    
    def test_run_emmernet(self):
        from locscale.emmernet.emmernet_functions import get_cubes, assemble_cubes
        from locscale.emmernet.run_emmernet import load_emmernet_model
        from locscale.include.emmer.ndimage.map_utils import load_map
        import numpy as np

        
        cubes_dictionary, cubes_array, filtered_signal_cubecenters = get_cubes(self.masked_emmap, self.stride, self.cube_size, self.mask)
        cube_1 = cubes_array[0]
        cube_size = cube_1.shape[0]
        batch_size = 8
    
        emmernet_model_1 = load_emmernet_model(self.inputs_dictionary)
        i=0
        cubes = cubes_array
        cubes_x = np.expand_dims(cubes, axis=4)
        cubes_predicted = np.empty((0, cube_size, cube_size, cube_size, 1))
        cubes_batch_X = np.empty((batch_size, cube_size, cube_size, cube_size, 1))
        cubes_batch_X = cubes_x[i:i+batch_size,:,:,:,:]

        ## Predict using model_based
        cubes_batch_predicted = emmernet_model_1.predict(x=cubes_batch_X, batch_size=batch_size, verbose=0)
        cubes_predicted = np.append(cubes_predicted, cubes_batch_predicted, axis=0)
        cubes_predicted = np.squeeze(cubes_predicted, axis=-1)

        self.assertTrue(cubes_predicted is not None)

        for predicted_cube in cubes_predicted:
            mean_predicted_cube = np.mean(predicted_cube)
            self.assertTrue(mean_predicted_cube < 15)
                        
        


    





        
        
            