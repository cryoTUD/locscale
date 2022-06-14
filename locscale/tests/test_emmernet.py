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

        
    def test_map_chunking(self):
        from locscale.emmernet.emmernet_functions import get_cubes, assemble_cubes
        import numpy as np
        print("Testing: map_chunking")

        emmap_shape = (252,252,252)
        cube_size = 32
        stride = 16
        emmap = np.random.randint(0,255,emmap_shape)

        # test 1: get_cubes
        cubes, cube_centers = get_cubes(emmap, stride, cube_size)

        # test 2: assemble_cubes
        
        assembled_cube = assemble_cubes(cubes, cube_centers, emmap_shape[0], cube_size)

        # test 3: check if assembled cube is same as original emmap
        self.assertTrue(np.allclose(assembled_cube, emmap))
        self.assertTrue(cubes[0].shape == (cube_size, cube_size, cube_size))
    
    def test_load_emmernet_model(self):
        from locscale.emmernet.run_emmernet import load_emmernet_model

        emmernet_type_1 = "model_based"
        emmernet_type_2 = "model_free"
        emmernet_type_3 = "ensemble"

        emmernet_model_1 = load_emmernet_model(emmernet_type_1, self.emmernet_model_folder)
        emmernet_model_2 = load_emmernet_model(emmernet_type_2, self.emmernet_model_folder)
        emmernet_model_3 = load_emmernet_model(emmernet_type_3, self.emmernet_model_folder)

        self.assertTrue(emmernet_model_1 is not None)
        self.assertTrue(emmernet_model_2 is not None)
        self.assertTrue(emmernet_model_3 is not None)
    
    def test_run_emmernet(self):
        from locscale.emmernet.emmernet_functions import get_cubes, assemble_cubes
        from locscale.emmernet.run_emmernet import load_emmernet_model
        from locscale.include.emmer.ndimage.map_utils import load_map
        import numpy as np

        emmap = np.random.normal(15,0.5,(252,252,252))
        cube_size = 32
        stride = 16
        batch_size = 8
        cubes, cube_centers = get_cubes(emmap, stride, cube_size)
        cube_1 = cubes[0]

        emmernet_model_1 = load_emmernet_model("model_based", self.emmernet_model_folder)
        i=0
        cubes = np.array(cubes)
        cubes_x = np.expand_dims(cubes, axis=4)
        print(cubes_x.shape)
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
                        
        


    





        
        
            