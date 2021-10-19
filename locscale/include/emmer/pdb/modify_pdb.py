#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 18 21:01:23 2021

@author: alok
"""
import numpy as np
import mrcfile
import gemmi
from locscale.include.emmer.ndimage.map_utils import convert_pdb_to_mrc_position

def set_pdb_cell_based_on_gradient(gemmi_st, emmap, apix, outpdb_name="modified.pdb", cell='occ'):
    gz, gx, gy = np.gradient(emmap)
    gradient_magnitude = np.sqrt(gx**2 + gy**2 + gz**2)
    gradient_magnitude_norm = gradient_magnitude / gradient_magnitude.max() * 100
    for model in gemmi_st:
        for chain in model:
            for res in chain:
                for atom in res:
                    pdb_pos = atom.pos.tolist()
                    mrc_pos = convert_pdb_to_mrc_position([pdb_pos], apix=apix)[0]
                    grad = gradient_magnitude_norm[mrc_pos[0], mrc_pos[1], mrc_pos[2]]
                    atom.occ = grad
    
    gemmi_st.write_pdb(outpdb_name)
                    