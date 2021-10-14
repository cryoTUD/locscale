#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep  9 11:34:06 2021

@author: alok
"""

import gemmi
import os
import numpy as np

folder = "/mnt/c/Users/abharadwaj1/Downloads/ForUbuntu/LocScale/tests/bfactor_correlation/emd3061/effect_of_add_blur/"


pseudo_atomic_model = folder+ "pdb5a63_shifted_refmac_refined.pdb"
atomic_model = folder+"pdb5a63_shifted_refmac_refined_blur10.pdb"

pseudo_gemmi_st = gemmi.read_structure(pseudo_atomic_model)

atomic_gemmi_st = gemmi.read_structure(atomic_model)

ns = gemmi.NeighborSearch(pseudo_gemmi_st[0], pseudo_gemmi_st.cell, 3).populate()  ## Search for corresponding atom in pseudo-atomic model

bfactor_correlation = {}
for model in atomic_gemmi_st:
    for chain in model:
        for res in chain:
            ca = res.get_ca()
            if ca is not None:
                bfactor_ca_atom = ca.b_iso
                
                atoms = ns.find_neighbors(ca, min_dist=0.1, max_dist=3)
                if len(atoms)>0:
                    dist = np.array([ca.pos.dist(atom.pos()) for atom in atoms])
                    index = dist.argmin()
                    nearest_atom = atoms[index]
                    chainID,resID,atID = nearest_atom.chain_idx, nearest_atom.residue_idx,nearest_atom.atom_idx
                    
                    pseudo_atom = pseudo_gemmi_st[0][chainID][resID][atID]
                    bfactor_pseudo_atom = pseudo_atom.b_iso
                    
                    
                    
                    if bfactor_pseudo_atom > 5:
                        bfactor_correlation[tuple(ca.pos.tolist())] = [bfactor_ca_atom, bfactor_pseudo_atom, dist[index]]

import pandas as pd
import seaborn as sns

df = pd.DataFrame(data=bfactor_correlation.values(), columns=['bf_ca','bf_pa','dist'])

df['diff'] = np.log(abs(df['bf_pa']-df['bf_ca']))

sns.relplot(data=df,x='bf_ca',y='bf_pa')

print("Correlation bf atomic vs pseudo-atomic: ",df['bf_ca'].corr(df['bf_pa'], method="spearman"))
    
                
            
            






