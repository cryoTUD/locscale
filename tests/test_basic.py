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
        from locscale.utils.file_tools import check_dependencies
        
        dependency=check_dependencies()      
        self.assertTrue(dependency)
