# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import mrcfile
import gemmi
import os
from locscale.include.emmer.ndimage.profile_tools import estimate_bfactor_standard, frequency_array, compute_radial_profile
from locscale.include.emmer.pdb.pdb_utils import get_coordinates
from locscale.include.emmer.ndimage.map_utils import convert_pdb_to_mrc_position
from locscale.include.emmer.ndimage.map_utils import extract_window
import random
import numpy as np

import seaborn as sns
import matplotlib.pyplot as plt


## User inputs
folder = "/mnt/c/Users/Alok/Downloads/ForUbuntu/tests/bfactor_correlation/emd5778"

emmap_path = os.path.join(folder, "using_atomic_model.mrc")
refined_model_path = os.path.join(folder, "pdb3j5p_refined_cropped.pdb")
local_window_size = 20  #A
sample_size = 500 # Number of different local windows to sample
fsc_resolution = 3.4

## Processing

emmap = mrcfile.open(emmap_path).data
apix = mrcfile.open(emmap_path).voxel_size.tolist()[0]

refined_model = gemmi.read_structure(refined_model_path)


## Find atomic locations

list_of_atoms = []
for model in refined_model:
    for chain in model:
        for res in chain:
            ca = res.get_ca()
            list_of_atoms.append(ca)
            
random_atoms = random.sample(list_of_atoms, sample_size)

# Find neighbors

ns = gemmi.NeighborSearch(refined_model[0], refined_model.cell, local_window_size).populate()
window_size_pixels = int(round(local_window_size/apix))

average_atomic_bfactor_model = []
average_bfactor_window = []

for random_atom in random_atoms:
    ## Find all neighbording atoms
    neighbors = ns.find_neighbors(random_atom, min_dist=0.1, max_dist=local_window_size)
    atoms = [refined_model[0][x.chain_idx][x.residue_idx][x.atom_idx] for x in neighbors]
    atomic_bfactor_list = np.array([x.b_iso for x in atoms])
    
    ## Find average bfactor based on average atomic bfactor
    average_atomic_bfactor = atomic_bfactor_list.mean()
    
    ## Find local bfactor based on local window
    
    # Find wilson_cutoff
    num_atoms = len(neighbors)
    print(num_atoms)
    mol_weight = num_atoms * 16  # daltons 
    wilson_cutoff_local = 1/(0.309 * np.power(mol_weight, -1/12))   ## From Amit Singer
    wilson_cutoff_local = np.clip(wilson_cutoff_local, apix*2*1.5, 10)
    ## Clipping necessary so that wilson_cutoff stays between reasonable values, can happen when there's too few atoms

    mrc_position = convert_pdb_to_mrc_position([random_atom.pos.tolist()], apix=apix)[0]
    local_window = extract_window(emmap, mrc_position, window_size_pixels)
    radial_profile = compute_radial_profile(local_window)
    freq = frequency_array(radial_profile, apix)
    bfactor_local_window = estimate_bfactor_standard(freq, radial_profile, wilson_cutoff=wilson_cutoff_local, fsc_cutoff=fsc_resolution)
    
    ## Compile results
    average_atomic_bfactor_model.append(average_atomic_bfactor)
    average_bfactor_window.append(abs(bfactor_local_window))
    

ax=sns.scatterplot(x=average_atomic_bfactor_model, y=average_bfactor_window)
ax.set(xlabel="Average atomic bfactors", ylabel="Bfactor estimation from slope of radial profile")
plt.savefig(os.path.join(folder, "avg_bfactor.jpg"))





    
    
    
    
    








