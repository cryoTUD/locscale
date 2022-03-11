#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar  9 17:45:00 2022

@author: alok
"""
#%% Introduction
'''
This script is used to find the effect of scattering atoms on the local bfactor correlation plots

Input: 
    1) refined  atomic model path
    2) number of scattering models (1A, 2A, 5A, 10A etc)
Output: 
    1) local bfactor correlation of atomic model map and scattered atomic model map
    2) local bfactor correlation of scattered atomic model map and atomic bfactors of scattered model
'''

import mrcfile
import gemmi
import os
import numpy as np
import random
from tqdm import tqdm
from locscale.include.emmer.pdb.pdb_tools import find_wilson_cutoff
from locscale.include.emmer.pdb.pdb_to_map import pdb2map
from locscale.include.emmer.pdb.pdb_utils import shake_pdb_within_mask
from locscale.include.emmer.pdb.pdb_tools import get_all_atomic_positions, get_atomic_bfactor_window
from locscale.include.emmer.ndimage.map_utils import extract_window, convert_pdb_to_mrc_position, resample_image, save_as_mrc
from locscale.include.emmer.ndimage.profile_tools import compute_radial_profile, frequency_array, estimate_bfactor_standard

import seaborn as sns
import matplotlib.pyplot as plt



#%% Inputs
folder = "/mnt/c/Users/abharadwaj1/Downloads/ForUbuntu/faraday_paper/Figures/get_perturbed_inputs"
refined_atomic_model_path = os.path.join(folder, "pdb6y5a_morphed_refined.pdb")
sample_mask_path = os.path.join(folder,"emd_10692_confidenceMap.mrc")
refined_model_gemmi_structure = gemmi.read_structure(refined_atomic_model_path)

rmsd_magnitudes = np.array([0.1, 0.5, 1, 2, 5, 10, 15,20])   ## in angstoerms


## Simulation parameters

apix = mrcfile.open(sample_mask_path).voxel_size.tolist()[0]   # angstoerm per pixel
map_shape = mrcfile.open(sample_mask_path).data.shape   ## (voxels, voxels, voxels)
boxsize_real_length = 25 ## angstoerm

boxsize = int(round(boxsize_real_length / apix))
nyquist_freq = apix*2

#%% Calculation

print("Calculating a reference model map")

refined_model_map = pdb2map(input_pdb=refined_model_gemmi_structure, apix=apix, size=map_shape, set_refmac_blur=True)

# Get several shaken gemmi structures
print("perturbed reference model with following magnitudes: {} ".format(rmsd_magnitudes))
shaken_structures = {}

for rmsd_magnitude in rmsd_magnitudes:
    shaken_structures[rmsd_magnitude] = shake_pdb_within_mask(refined_atomic_model_path, sample_mask_path, rmsd_magnitude=rmsd_magnitude,use_pdb_mask=True)
    shaken_structures[rmsd_magnitude].write_pdb(os.path.join(folder, "pdb6y5a_rmsd_{}_pm_perturbed_using_mask.pdb".format(int(rmsd_magnitude*100))))



# Basic calculation for estimating bfactors
gemmi_model = refined_model_gemmi_structure[0]
num_atoms = gemmi_model.count_atom_sites()
wilson_cutoff = find_wilson_cutoff(num_atoms=num_atoms)

shaken_maps = {}

print("Calculating simulated maps from the shaken models")
for rmsd_magnitude in rmsd_magnitudes:
    shaken_maps[rmsd_magnitude] = pdb2map(input_pdb=shaken_structures[rmsd_magnitude], apix=apix, size=map_shape, set_refmac_blur=True)

for rmsd_magnitude in rmsd_magnitudes:
    save_as_mrc(shaken_maps[rmsd_magnitude], output_filename=os.path.join(folder, "perturbed_map_rmsd_{}_pm.mrc".format(int(rmsd_magnitude*100))), apix=apix)

#%%
# Find atomic positions from main structure

atomic_positions_main = get_all_atomic_positions(refined_model_gemmi_structure)
num_positions = min(1000, 0.5 * num_atoms)
random_sample = random.sample(list(atomic_positions_main), num_positions)

random_centers_real_unit = []
for center_real in random_sample:
    center_r = tuple([center_real[0]+np.random.uniform(-3,3),center_real[1]+np.random.randint(-3,3),center_real[2]+np.random.randint(-3,3)])
    random_centers_real_unit.append(center_r)
    
# Find radial profiles 

bfactor_using_maps = {}
bfactor_using_model = {}
quality_of_fit_analysis = {}

center_int = []
for center_real in random_centers_real_unit:
    center = tuple(convert_pdb_to_mrc_position([center_real], apix)[0])
    center_int.append(center)
    
print("Calculating bfactors")
for center in tqdm(center_int, desc="Calculating local bfactors from maps"):
    
    temp = {}
    main_window = extract_window(refined_model_map, center, boxsize)
    ref_radial_profile = compute_radial_profile(main_window)
    freq = frequency_array(ref_radial_profile, apix=apix)
    
    local_wilson_cutoff = find_wilson_cutoff(num_atoms=ref_radial_profile[0])
    quality_of_fit = {}
    bfactor_main_window, r2_main = estimate_bfactor_standard(freq, ref_radial_profile, wilson_cutoff=local_wilson_cutoff, fsc_cutoff=nyquist_freq, return_fit_quality=True, standard_notation=True)
    shaken_windows = {}
    shaken_radial_profile = {}
    shaken_window_bfactor = {}
    shaken_window_bfactor[0] = bfactor_main_window
    quality_of_fit[0] = r2_main
    for rmsd_magnitude in rmsd_magnitudes:
        shaken_windows[rmsd_magnitude] = extract_window(shaken_maps[rmsd_magnitude], center, boxsize)
        shaken_radial_profile[rmsd_magnitude] = compute_radial_profile(shaken_windows[rmsd_magnitude])
        shaken_window_bfactor[rmsd_magnitude], quality_of_fit[rmsd_magnitude] = estimate_bfactor_standard(freq, shaken_radial_profile[rmsd_magnitude], 
                                                                             wilson_cutoff=local_wilson_cutoff, fsc_cutoff=nyquist_freq, return_fit_quality=True, standard_notation=True)
        
    
    bfactor_using_maps[center] = shaken_window_bfactor
    quality_of_fit_analysis[center] = quality_of_fit
#%%
for atomic_position in tqdm(random_centers_real_unit, desc="Calculating local bfactor from model"):
    shaken_bfactor_model = {}
    
    shaken_bfactor_model[0] = get_atomic_bfactor_window(refined_model_gemmi_structure, atomic_position, boxsize_real_length)
    
    for rmsd_magnitude in rmsd_magnitudes:
        shaken_bfactor_model[rmsd_magnitude] = get_atomic_bfactor_window(shaken_structures[rmsd_magnitude], atomic_position, boxsize_real_length)
    
    bfactor_using_model[tuple(atomic_position)] = shaken_bfactor_model

import pandas as pd

df = pd.DataFrame()

for center in center_int:
    
    rmsd = [0]
    for r in rmsd_magnitude:
        rmsd.append(r)
    
    temp_list = []
    for r in rmsd:
        temp_list.append(bfactor_using_maps[center][r])
    
    for r in rmsd:
        temp_list.append(bfactor_using_model[tuple(center_real)][r])
    
    for r in rmsd:
        temp_list.append(quality_of_fit_analysis[center][r])

    df[center] = temp_list

columns = []
for r in rmsd:
    columns.append("shaken_map_{}_A".format(r))

for r in rmsd:
    columns.append("shaken_model_{}_A".format(r))
    
for r in rmsd:
    columns.append("quality_fit_{}_A".format(r))

df.index = columns
df = df.T

plt.figure(1)
for i in [0]+rmsd_magnitude:
    avg_qfit = df['quality_fit_{}_A'.format(i)].mean()
    plt.plot(i, avg_qfit, 'ko')
    plt.xlabel("RMSD (A)")
    plt.ylabel("Average quality of fit")
    

plt.figure(2)

correlations = []
for i in rmsd_magnitude:
    correlation = df['shaken_map_{}_A'.format(i)].corr(df["shaken_map_0_A"])
    #plt.plot(i, abs(correlation), 'ko')
    #plt.xlabel("RMSD (A)")
    #plt.ylabel("Correlation with unshaken atomic model map")
    correlations.append(correlation)
    
    
#%%

sns.set_theme(context="paper", font="Helvetica", font_scale=1.5)
sns.set_style("white")
kwargs = dict( marker='o', alpha=0.3)

sns.relplot(data=df, x="shaken_map_0_A",y="shaken_map_40_A", **kwargs)
plt.xlabel("Guinier B-factor unperturbed model map ($\AA^{-2}$)")
plt.ylabel("Guinier B-factor for RMSD: 10 $\AA$  ($\AA^{-2}$)")
plt.tight_layout()
#plt.savefig("Figure_2_bfactor.svg", dpi=600)


        
#%%

sns.lineplot(x=rmsd_magnitude, y=correlations, marker="o", linewidth=3, markersize=12)
plt.xlabel("RMSD ($\AA$)")
plt.ylabel("Correlation Coefficient")
plt.tight_layout()
#plt.savefig("Correlation_coeffient.svg", dpi=600)

#%%

# Save maps and pdbs
for rmsd in rmsd_magnitude:
    st = shaken_structures[rmsd]
    simmap = shaken_maps[rmsd]
    st.write_pdb(os.path.join(folder, "pdb6y5a_rmsd_{}_A.pdb".format(rmsd)))
    
    save_as_mrc(simmap, output_filename=os.path.join(folder, "perturbed_map_rmsd_{}_A.mrc".format(rmsd)), apix=apix)

    
    