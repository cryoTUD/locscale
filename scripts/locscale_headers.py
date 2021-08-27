#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 17 11:31:51 2020

@author: alok
"""
import argparse
from pam_headers import *
import ipykernel
from pseudomodel_analysis import *
import pandas as pd
from emmer.pdb.pdb_tools import *
import pprint

paths = os.environ['PATH']
allpaths = paths.split(':')
for path in allpaths:
    if 'ccpem' in path and 'bin' in path:
        path_to_ccpem = path[:-4]
    if 'ccp4' in path and 'bin' in path:
        path_to_ccp4 = path[:-4]
    if 'locscale' in path and 'scripts' in path:
        path_to_locscale = path[:-8]

path_to_ccpem = "/home/alok/soft/ccpem-20200424"
path_to_ccp4="/home/alok/soft/ccp4-7.1/ccp4-7.1"

def run_FDR(emmap_path,window_size,fdr=0.01,verbose=True,filter_cutoff=None):
    '''
    

    Parameters
    ----------
    emmap_path : string
        Path to the EM Map which needs thresholding. Example: 'path/to/map.mrc'
    window_size : int
        Window size required for FDR thresholding
    verbose : bool, optional
        Print statistics if True. The default is True.

    Returns
    -------
    mask_path : string
        path to mask file. Example: 'path/to/mask.mrc'
    mask_mrc : mrfile.mrc() 
        mrcfile object which contains volume and header information about the mask

    '''
    import os, sys
    
    # Preprocessing EM Map Path
    
    # Apply filter if filter_cutoff is not None

    if verbose:
        print("Now starting FDR procedure using the following parameters: \n"
                 "Window size: "+str(window_size)+"\n"
                 "Filter cutoff: \n"+str(filter_cutoff))
    
                    
    try:
        ## First attempt to use Emmer package to compute FDR map
        from emmer.ndimage.map_utils import average_voxel_size, compute_FDR_confidenceMap_easy, save_as_mrc
        emmap = mrcfile.open(emmap_path).data
        voxel_size_record = mrcfile.open(emmap_path).voxel_size
        apix = average_voxel_size(voxel_size_record)
        fdr = fdr
        fdr_mask, fdr_threshold = compute_FDR_confidenceMap_easy(
            emmap, apix=apix, fdr=fdr, window_size=window_size, 
            lowPassFilter_resolution=filter_cutoff)
        
        print("FDR threshold found to be: \t", fdr_threshold)
        emmap_path_without_ext = emmap_path[:-4]
        mask_path = emmap_path_without_ext + "_confidenceMap.mrc"
        
        save_as_mrc(fdr_mask, output_filename=mask_path, 
                    apix=voxel_size_record.tolist(), origin=0)
        
        if os.path.exists(mask_path):
            if verbose:
                print("FDR Procedure completed. \n"+
                      "Mask path: "+mask_path+"\n")
            
            return mask_path
        else:
            print("FDR process failed. Returning none")
            return None
        
    
    except:    
        print(sys.exc_info())
        print("Could not use the FDRUtil python package. Reverting to ccpem version of FDRUtil")
        path_to_FDR_script = path_to_ccpem+"/lib/py2/FDRcontrol.pyc"
        fdr_command_line = "ccpem-python "+path_to_FDR_script+" --em_map "+emmap_path+" -method BY --testProc rightSided --window_size "+str(window_size)
                   
        fdr_output = run(fdr_command_line.split(),stdout=PIPE)
        emmap_name = emmap_path[:-4]
        mask_path = emmap_name+"_confidenceMap.mrc"
        if os.path.exists(mask_path):
            mask_mrc = mrcfile.open(mask_path)
            if verbose:
                print("FDR Procedure successful! \n"+
                      "Mask Path: "+mask_path)
            return mask_path 
        else:
            print("FDR procedure unsuccessful. Returning None")
            return None
            

def measure_mask_parameters(mask_path,protein_density=1.35,average_atomic_weight=13.14,verbose=True,detailed_report=False):
    '''
    Parameters
    ----------
    mask_path : string 
        Path to mask file
    edge_threshold : float 
        The threshold to strictly binarize the FDR map at the edges
    protein_density : float, optional
        Average protein density to calculate number of atoms. The default is 1.35.
    average_atomic_weight : float, optional
        Atomic weight of an "average atom present in protein". 
        The default is 13.14. Found using calculating average atomic weight (molecular mass/num atoms) for 500 proteins.
    verbose : bool, optional
        Print statistics if True. The default is True.

    Returns
    -------
    num_atoms : int
        Estimated number of atoms based on mask volume, protein density and average atomic weight

    '''
    from scipy.constants import Avogadro
    from emmer.ndimage.map_utils import average_voxel_size
    
    mask_mrc = mrcfile.open(mask_path)
    voxelsize = average_voxel_size(mask_mrc.voxel_size)
    
    ang_to_cm = 1e-8
    mask = mask_mrc.data
    mask.setflags(write=1)
    mask[mask<0.5] = 0
    mask[mask>=0.5] = 1
    mask_vol = mask.sum()*(voxelsize*ang_to_cm)**3
    mask_vol_A3 = mask.sum()*voxelsize**3
    
    protein_mass = protein_density * mask_vol
    num_moles = protein_mass / average_atomic_weight
    num_atoms = int((num_moles * Avogadro).round())
    maskshape = mask.shape
    mask_dims = [maskshape[0]*voxelsize,maskshape[1]*voxelsize,maskshape[2]*voxelsize]
        
    if verbose:
        print("Mask parameters calculated are: \n"+
              "Mask volume: "+str(round(mask_vol_A3,3))+" A$^3$ \n"+
              "Protein mass: "+str(round(1e21*protein_mass))+" zg\n"+
              "Num atoms: "+str(num_atoms)+"\n") 
        
    if not detailed_report:
        return num_atoms,mask_dims
    else:
        return mask_vol_A3, protein_mass, num_atoms, mask_dims,maskshape

def create_pseudo_mask(emmap_path,threshold):
    emmap_mrc = mrcfile.open(emmap_path)
    emmap = emmap_mrc.data
    emmap.setflags(write=1)
    emmap[emmap<threshold] = 0
    emmap[emmap>=threshold] = 1
    mask_path = emmap_path[:-4]+'_mask.mrc'
    with mrcfile.new(mask_path,overwrite=True) as maskmrc:
        maskmrc.set_data(emmap)
        maskmrc.voxel_size = emmap_mrc.voxel_size
        maskmrc.header.origin = emmap_mrc.header.origin
    return mask_path

def run_pam(emmap_path,mask_path,threshold,num_atoms,method,bl,
            g=None,friction=None,scale_map=None,scale_lj=None,total_iterations=100,verbose=True):
    '''
    

    Parameters
    ----------
    emmap_path : string
        Path to the EM Map which contains the pseudo-atomic model
    mask_path : string
        Path to the Masked map 
    threshold : float
        Threshold required to strictly binarize the FDR masked map, especially at the edges
    num_atoms : int
        Number of atoms to fill in the pseudo-atomic model
    method : string
        Method to generatre the pseudo-atomic model. Value is either: 
            'gradient' for high resolution EM Maps
            'random_placement_with_kick' for poor resolution maps 
    bl : float
        bl = "bond length" refers to the minimum distance between two atoms in the pseudo-atomic model which needs to be satisfied
    
    total_iterations : int
        Total number of iterations to run
        
    --- following parameters required if method = 'gradient'    
    
    g : float
        Scale acceleration due to gradient forces
    
    friction : float
        Friction coefficient for solver to converge
    
    scale_map : float
        For overall scaling of the gradient potential. Default is 1.
    scale_lj : float
        For overall scaling of the inter-atomic forces, modelled by LJ Potential. Default is 1. 
        
    

    Returns
    -------
    pseudomodel_path : string
        Path of the output pseudo-atomic model 
        

    '''
    
    mrc = mrcfile.open(emmap_path)
    emmap = mrc.data
    voxelsize = mrc.voxel_size.x

    mask = mrcfile.open(mask_path).data
    pseudomodel = extract_model_from_mask(mask,num_atoms,threshold=threshold)
    
    emmap_shape = emmap.shape
    unitcell = gemmi.UnitCell(emmap_shape[0]*voxelsize,emmap_shape[1]*voxelsize,emmap_shape[2]*voxelsize,90,90,90)
    
    if verbose:
        print("Running pseudoatomic model generator to add "+str(num_atoms)+" atoms inside the volume using the method: "+method)
    if method=='gradient':
        if g is None:
            g = 10
        if scale_lj is None:
            scale_lj = 1
        if scale_map is None:
            scale_map = 1
        if friction is None:
            friction = 10
            
        gz,gy,gx = np.gradient(emmap)
        arranged_points = main_solver3D(
            emmap,gx,gy,gz,pseudomodel,g=g,friction=friction,min_dist_in_angst=bl,voxelsize=voxelsize,dt=0.1,capmagnitude_lj=100,epsilon=1,scale_lj=scale_lj,
            capmagnitude_map=100,scale_map=scale_map,total_iterations=total_iterations, path_for_gemmi_models=None,emmap_path=None,mask_path=None,returnPointsOnly=True,
            integration='verlet',verbose=True)
        mask_name = mask_path[:-4]
        pseudomodel_path = mask_name+"_gradient_pseudomodel.pdb"

    elif method=='random' or method=='kick' or method == 'random_placement_with_kick':
        arranged_points = main_solver_kick(
                pseudomodel,min_dist_in_angst=bl,voxelsize=voxelsize,total_iterations=99,returnPointsOnly=True,verbose=True)
        mask_name = mask_path[:-4]
        pseudomodel_path = mask_name+"_kick_pseudomodel.pdb"
    
    arranged_points.write_pdb(pseudomodel_path,voxelsize=voxelsize,unitcell=unitcell)
    
    
    if os.path.exists(pseudomodel_path):    
        print("The location of the pseudomodel generated is: "+pseudomodel_path+'\n\n')
        return pseudomodel_path
    else:
        print("uhhu, something wrong with the pseudomodel generator! Returning None")        
        return None
    
        

def run_refmac(model_path,model_name,map_path,resolution,maskdims,verbose=True):
    
    path_to_run_refmac = path_to_locscale+"/scripts/run_refmac.zsh"
    refmac_command_line = "zsh "+path_to_run_refmac+" "+model_path+" "+model_name+" "+map_path+" "+str(round(resolution,2))+" "+path_to_ccpem+" "+path_to_ccp4+" "+str(maskdims[0])+" "+str(maskdims[1])+" "+str(maskdims[2])
    if verbose:
        print("Running REFMAC to refine the pseudo-atomic model using \n"+
              "Path to run_refmac: "+path_to_run_refmac+
              "Command line: "+refmac_command_line)
        
    refmac_output = run(refmac_command_line.split(),stdout=PIPE)
    refined_model_path = model_name+"_refmac_refined.pdb"
        
    if os.path.exists(refined_model_path):
        if verbose: 
            print("The refined PDB model is: "+refined_model_path+"\n\n")    
        return refined_model_path
    else:
        print("Uhhoh, something wrong with the REFMAC procedure. Returning None")
        return None
    
    
    
def run_refmap(model_path,emmap_path,mask_path,resolution=None,verbose=True):
    '''
    Function to obtain reference map using structure factors determined by atomic model.
    This function uses gemmi.DensityCalculatorE() function to calculate a grid of intensities.
    
    Required modules: emmer

    Parameters
    ----------
    model_path : str
        path for atomic model 
    emmap_path : str
        path to emmap map
    mask_path : str
        path to mask 
    verbose : bool, optional
        The default is True.

    Returns
    -------
    reference_map : str
        Path to reference map generated.

    '''
    from emmer.pdb.pdb_to_map import pdb2map
    from emmer.ndimage.map_utils import average_voxel_size, save_as_mrc, read_gemmi_map, compare_gemmi_grids
    from emmer.ndimage.map_tools import get_center_of_mass
    import pandas as pd
    
    if verbose: 
        print("Now simulating Reference Map using Refined Atomic Model")
    
    # Read inputs from filesystem
    
    emmap_data, grid_input = read_gemmi_map(emmap_path, return_grid=True)
    mask = read_gemmi_map(mask_path)
    pdb_structure = gemmi.read_structure(model_path)
    
    
    ## Generate parameters of the simulated map
    voxelsize = grid_input.spacing   
    unitcell = grid_input.unit_cell
    
    ## Simulate a reference map from the input atomic model in the pdb_structure variable
    
    refmap_data, grid_simulated = pdb2map(input_pdb=pdb_structure, unitcell=unitcell, size=emmap_data.shape,
                                          return_grid=True, align_output=True)
    
    ## Output filename
    reference_map_path = model_path[:-4]+"_4locscale.mrc"
    save_as_mrc(map_data=refmap_data,output_filename=reference_map_path, apix=grid_simulated.spacing, origin=0)   
    
    ## Checklist: 
    
    correlation = compute_real_space_correlation(mrcfile.open(emmap_path).data, 
                                                 mrcfile.open(reference_map_path).data)
    
    grid_comparison = compare_gemmi_grids(read_gemmi_map(emmap_path, return_grid=True)[1], 
                                          read_gemmi_map(reference_map_path,return_grid=True)[1])
    # Since grid comparison and correlation are critical, they are done on the saved filesystem 
    # and not on the files in memory. This is to avoid any errors that might have happened during 
    # save_as_mrc operation
    
    center_of_mass_experimental = get_center_of_mass(emmap_data,apix=grid_input.spacing)
    center_of_mass_simulated = get_center_of_mass(refmap_data,apix=grid_simulated.spacing)
    center_of_mass_atomic_model = pdb_structure[0].calculate_center_of_mass().tolist()
    
    reporter = {}
    reporter['Model-map_Correlation'] = correlation    
    reporter['COM:_Experimental_map'] = center_of_mass_experimental    
    reporter['COM:_Simulated_map'] = center_of_mass_simulated    
    reporter['COM:_Atomic_model'] = center_of_mass_atomic_model
    reporter['Grid comparison'] = grid_comparison['final'].all()
    
    ## Add checkpoint: center of mass of pseudo-model, simulated map and original map, (2) correlation (3) Axis order
    if os.path.exists(reference_map_path):
        if verbose: 
            
            print("The reference map is at: "+reference_map_path+"\n\n")
            pprint.pprint(reporter)
            
        return reference_map_path
    else:
        print("Reference map was not generated. Returning none")
        return None
    