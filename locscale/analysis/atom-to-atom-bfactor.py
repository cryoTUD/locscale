#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep  9 11:34:06 2021

@author: alok
"""

from myheaders import *
folder = "/mnt/c/Users/abharadwaj1/Downloads/ForUbuntu/LocScale/temp_files/"

emmap,grid = read_gemmi_map(folder+"emd5778_unfiltered.mrc", return_grid=True)

pseudo_atomic_model = folder+ "3j5p_pseudomodel_gradient_refmac_refined_iter45.pdb"
atomic_model = folder+"pdb3j5p_refined_cropped.pdb"

pseudo_gemmi_st = gemmi.read_structure(pseudo_atomic_model)

atomic_gemmi_st = gemmi.read_structure(atomic_model)

ns = gemmi.NeighborSearch(pseudo_gemmi_st[0], pseudo_gemmi_st.cell, 3).populate()  ## Search for corresponding atom in pseudo-atomic model

bfactor_correlation = {}
for model in atomic_gemmi_st:
    for chain in model:
        for res in chain:
            ca = res.get_ca()
            bfactor_ca_atom = ca.b_iso
            
            atoms = ns.find_neighbors(ca, min_dist=0.1, max_dist=3)
            if len(atoms)>0:
                dist = np.array([ca.pos.dist(atom.pos()) for atom in atoms])
                index = dist.argmin()
                nearest_atom = atoms[index]
                chainID,resID,atID = nearest_atom.chain_idx, nearest_atom.residue_idx,nearest_atom.atom_idx
                
                pseudo_atom = pseudo_gemmi_st[0][chainID][resID][atID]
                bfactor_pseudo_atom = pseudo_atom.b_iso
                
                map_value = grid.interpolate_value(ca.pos)
                
                if bfactor_pseudo_atom > 5:
                    bfactor_correlation[tuple(ca.pos.tolist())] = [bfactor_ca_atom, bfactor_pseudo_atom, dist[index], map_value]

df = pd.DataFrame(data=bfactor_correlation.values(), columns=['bf_ca','bf_pa','dist','map_value'])

df['diff'] = np.log(abs(df['bf_pa']-df['bf_ca']))

sns.relplot(data=df,x='bf_ca',y='bf_pa',hue='map_value')

print("Correlation bf atomic vs pseudo-atomic: ",df['bf_ca'].corr(df['bf_pa'], method="spearman"))
    
                
            
            






