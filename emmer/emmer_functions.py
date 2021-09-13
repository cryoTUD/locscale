#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 14 19:56:16 2021

@author: alok
"""

import emmer
from inspect import getmembers, isfunction, ismodule
from pprint import pprint

def get_module_names(module_object):
    modulenames = {}
    module_members = getmembers(module_object, ismodule)
    for module in module_members:
        modulenames[module[0]] = module[0]
    
    return modulenames

def getfunctions(module_object):
    functions = {}
    function_members = getmembers(module_object, isfunction)
    for function in function_members:
        functions[function[0]] = function[0]
    
    return functions
    
def getmodules(module_object):
    modules = {}
    module_members = getmembers(module_object, ismodule)
    
    for module in module_members:
        directory = {}
        if module[0] not in ['headers']:
            for submodule in getmembers(module[1],ismodule):
                submodule_dictionary = getfunctions(submodule[1])
                directory[submodule[0]] = submodule_dictionary
            modules[module[0]] = directory
        else:
            modules[module[0]] = get_module_names(module[1])
    
    return modules

emmer_module = getmodules(emmer)


#outfile = open("modules_list_full.txt","w")
#pprint(emmer_module, stream=outfile)
#outfile.close()

pprint(emmer_module)
