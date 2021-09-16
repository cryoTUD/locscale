# -*- coding: utf-8 -*-
"""
Created on Thu Jul 29 10:16:11 2021

"""

import argparse
msg = "loads in model and target map and sets the unit cell of the model to that of the target map."
parser = argparse.ArgumentParser(prog="set unit cell of pdb to target map", description=msg)
parser.add_argument("Model", help="input model. Must be .pdb file")
parser.add_argument("Map", help="target map to set model unit cell to")
parser.add_argument("-o","--output", help="path to and name of new pdb file with updated unit cell")
args = parser.parse_args()
from emmer.pdb.pdb_tools import add_cryst1_line

add_cryst1_line(pdb_path=args.Model,
                emmap_path=args.Map,
                new_pdb_path=args.output)



