#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Sep 12 18:55:43 2021

@author: alok
"""

## Unit tests to check scaling 
import unittest

class TestBasic(unittest.TestCase):
    
    def test_dependencies(self):
        from locscale.pseudomodel.pseudomodel_headers import check_dependencies
        
        dependency=check_dependencies()
        print(dependency)
        keys = list(dependency.keys())
        check = False
        if 'ccpem' in keys and 'ccp4'in keys and 'locscale' in keys:
            check = True
        self.assertTrue(check)
