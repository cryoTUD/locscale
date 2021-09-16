#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 17 11:31:51 2020

@author: alok
"""

import numpy as np


def check_dependencies():
    import os
    dependency = {}
    paths = os.environ['PATH']
    allpaths = paths.split(':')
    
    for path in allpaths:
        if 'ccpem' in path and 'bin' in path:
            path_to_ccpem = "/".join(path.split("/")[:-1])
            dependency['ccpem'] = path_to_ccpem
        if 'ccp4' in path and 'bin' in path:            
            path_to_ccp4 = "/".join(path.split("/")[:-1])
            dependency['ccp4'] = path_to_ccp4
        if 'locscale' in path.split("/"):
            locscale_index=path.split("/").index("locscale")
            path_to_locscale = "/".join(path.split("/")[:locscale_index+1])
            dependency['locscale'] = path_to_locscale
    
    '''
    if path_to_ccp4 is None or path_to_ccpem is None or path_to_locscale is None:
        print("Required dependencies not found! ")
        print("CCPEM Location at: ",path_to_ccpem)
        print("CCP4 Location at: ",path_to_ccp4)
        print("LocScale Location at: ",path_to_locscale)
        
        exit
    else:
        print("Required dependencies are found at: \n ")
        print("CCPEM Location at: ",path_to_ccpem)
        print("CCP4 Location at: ",path_to_ccp4)
        print("LocScale Location at: ",path_to_locscale)
    '''  
    
    return dependency
    
def number_of_segments(fsc_resolution):
    if fsc_resolution < 3:
        return 4
    elif fsc_resolution >= 3:
        return 3
    elif fsc_resolution >= 5:
        return 2
    else:
        print("Warning: resolution too low to estimate cutoffs. Returning 1")
        return 1
    

def prepare_sharpen_map(emmap_path,wilson_cutoff,fsc_resolution,add_blur=0,return_processed_files=False, output_file_path=None):
    from emmer.ndimage.profile_tools import compute_radial_profile, estimate_bfactor_through_pwlf, frequency_array
    from emmer.ndimage.map_utils import average_voxel_size, save_as_mrc
    from emmer.ndimage.map_tools import sharpen_maps, estimate_global_bfactor_map
    import mrcfile
    
    emmap_mrc = mrcfile.open(emmap_path)
    emmap_unsharpened = emmap_mrc.data
    apix=average_voxel_size(emmap_mrc.voxel_size)
    
    bfactor, z, slopes, fit = estimate_global_bfactor_map(emmap=emmap_unsharpened, apix=apix, wilson_cutoff=wilson_cutoff, fsc_cutoff=fsc_resolution)
    
    print("bfactor: {:.3f}, breakpoints: {} and slopes: {}".format(bfactor, (1/np.sqrt(z)).round(2),slopes))
    if add_blur != 0:
        bfactor  += add_blur  ## Use add_blur if you wanna add blur to the emmap before refining
        print("Final overall bfactor of emmap expected to be {:.2f}".format(-1*add_blur))
        
    sharpened_map = sharpen_maps(emmap_unsharpened, apix=apix, global_bfactor=bfactor)
    
    bfactor_final, z_final, slopes_final, _ = estimate_global_bfactor_map(emmap=sharpened_map, apix=apix, wilson_cutoff=wilson_cutoff, fsc_cutoff=fsc_resolution)
    
    print("Final overall bfactor of emmap computed to be {:.2f}".format(bfactor_final))
    print("bfactor: {:.3f}, breakpoints: {} and slopes: {}".format(bfactor_final, (1/np.sqrt(z_final)).round(2),slopes_final))
        
    
    output_filename = emmap_path[:-4] +"_global_sharpened.mrc"
    save_as_mrc(map_data=sharpened_map, output_filename=output_filename, apix=apix, origin=0)
    
    if return_processed_files:
        print("Returning: sharpend_map_path, [rp_unsharp, rp_sharp, bfactor]")
        return output_filename, fit
    else:
        return output_filename

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
    from subprocess import run, PIPE
    import mrcfile
    # Preprocessing EM Map Path
    
    # Apply filter if filter_cutoff is not None

    if verbose:
        print("Now starting FDR procedure using the following parameters: \n"
                 "Window size: "+str(window_size)+"\n"
                 "Filter cutoff: "+str(filter_cutoff))
    
                    
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
        path_to_ccpem = check_dependencies()['ccpem']
        path_to_FDR_script = path_to_ccpem+"/lib/py2/FDRcontrol.pyc"
        fdr_command_line = "ccpem-python "+path_to_FDR_script+" --em_map "+emmap_path+" -method BY --testProc rightSided --window_size "+str(window_size)
                   
        run(fdr_command_line.split(),stdout=PIPE)
        emmap_name = emmap_path[:-4]
        mask_path = emmap_name+"_confidenceMap.mrc"
        if os.path.exists(mask_path):
            if verbose:
                print("FDR Procedure successful! \n"+
                      "Mask Path: "+mask_path)
            return mask_path 
        else:
            print("FDR procedure unsuccessful. Returning None")
            return None
            
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
    import os
    import mrcfile
    import gemmi
    from locscale.pseudomodel.pseudomodel_solvers import main_solver3D, main_solver_kick
    from locscale.pseudomodel.pseudomodel_classes import extract_model_from_mask
    
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
        gz,gy,gx = np.gradient(emmap)
        masked_grad_magnitude = mask * np.sqrt(gx**2 + gy**2 + gz**2)
        max_gradient = masked_grad_magnitude.max()
        if g is None:
            g = round(100 / max_gradient)
        if scale_lj is None:
            scale_lj = 1
        if scale_map is None:
            scale_map = 1
        if friction is None:
            friction = 10
            
        
        arranged_points = main_solver3D(
            emmap,gx,gy,gz,pseudomodel,g=g,friction=friction,min_dist_in_angst=bl,voxelsize=voxelsize,dt=0.1,capmagnitude_lj=100,epsilon=1,scale_lj=scale_lj,
            capmagnitude_map=100,scale_map=scale_map,total_iterations=total_iterations, compute_map=None,emmap_path=None,mask_path=None,returnPointsOnly=True,
            integration='verlet',verbose=verbose)
        mask_name = mask_path[:-4]
        pseudomodel_path = mask_name+"_gradient_pseudomodel.pdb"

    elif method=='random' or method=='kick' or method == 'random_placement_with_kick':
        arranged_points = main_solver_kick(
                pseudomodel,min_dist_in_angst=bl,voxelsize=voxelsize,total_iterations=total_iterations,returnPointsOnly=True,verbose=verbose)
        mask_name = mask_path[:-4]
        pseudomodel_path = mask_name+"_kick_pseudomodel.pdb"
    
    arranged_points.write_pdb(pseudomodel_path,voxelsize=voxelsize,unitcell=unitcell)
    
    
    if os.path.exists(pseudomodel_path):    
        print("The location of the pseudomodel generated is: "+pseudomodel_path+'\n\n')
        return pseudomodel_path
    else:
        print("uhhu, something wrong with the pseudomodel generator! Returning None")        
        return None
    
def is_pseudomodel(input_pdb_path):
    '''
    Function to check if a pdb at a pdb_path is a pseudo_atomic model based on number of water molecules

    Parameters
    ----------
    input_pdb_path : TYPE
        DESCRIPTION.

    Returns
    -------
    is_pseudomodel_check : bool

    '''
    import gemmi
    gemmi_st = gemmi.read_structure(input_pdb_path)
    
    num_atoms = gemmi_st[0].count_atom_sites()
    
    num_waters = 0
    for model in gemmi_st:
        for chain in model:
            for res in chain:
                for atom in res:
                    if atom.element.name == 'O':
                        num_waters += 1
    
    if num_waters == num_atoms:
        print("Number of waters {} equal to number of atoms {}. \nHence {} is a pseudo-atomic model".format(num_waters, num_atoms, input_pdb_path))
        return True
        print("Number of waters {} not equal to number of atoms {}.\nHence {} is a pseudo-atomic model".format(num_waters, num_atoms, input_pdb_path))
        return False
    


def run_refmac(model_path,map_path,resolution,  num_iter,only_bfactor_refinement,verbose=True):
    import os
    from subprocess import run, PIPE
    from emmer.pdb.pdb_utils import get_bfactors
    import mrcfile
    
    path_to_locscale = check_dependencies()['locscale']
    path_to_ccpem = check_dependencies()['ccpem']
    path_to_ccp4 = check_dependencies()['ccp4']
    
    if only_bfactor_refinement:
        path_to_run_refmac = path_to_locscale+"/locscale/utils/run_refmac.zsh"
    else:
        path_to_run_refmac = path_to_locscale+"/locscale/utils/run_refmac_restrained.zsh"
        
    
    model_name = model_path[:-4]
    emmap_mrc = mrcfile.open(map_path)
    map_dims = emmap_mrc.header.cella.tolist()
    
    refmac_command_line = "zsh "+path_to_run_refmac+" "+model_path+" "+model_name+" "+map_path+" "+str(round(resolution,2))+" "+path_to_ccpem+" "+path_to_ccp4+" "+str(map_dims[0])+" "+str(map_dims[1])+" "+str(map_dims[2])+" "+str(num_iter)
    if verbose:
        print("Running REFMAC to refine the pseudo-atomic model using \n"+
              "Path to run_refmac: "+path_to_run_refmac+"\n"+
              "Command line: \n"+refmac_command_line)
        
    refmac_output = run(refmac_command_line.split(),stdout=PIPE)
    refined_model_path = model_name+"_refmac_refined.pdb"
    bfactors = get_bfactors(in_model_path=refined_model_path)
    
    if os.path.exists(refined_model_path):
        if verbose: 
            print("The refined PDB model is: "+refined_model_path+"\n\n")    
            print("B factor range: \t ({:.2f} to {:.2f})".format(min(bfactors),max(bfactors)))
            
        return refined_model_path
    else:
        print("Uhhoh, something wrong with the REFMAC procedure. Returning None")
        return None
    
def normalise_intensity_levels(from_emmap, to_levels=[0,1]):
    normalise_between_zero_one = (from_emmap - from_emmap.min()) / (from_emmap.max() - from_emmap.min())
    to_levels = np.array(to_levels)
    
    min_value = to_levels.min()
    max_value = to_levels.max()
    scale_factor = max_value-min_value
    
    normalised = min_value + normalise_between_zero_one * scale_factor
    
    return normalised    
    
def run_refmap(model_path,emmap_path,mask_path,add_blur=0,resolution=None,verbose=True):
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
    import os
    import gemmi
    import mrcfile
    import pprint
    from emmer.pdb.pdb_to_map import pdb2map
    from emmer.ndimage.map_utils import average_voxel_size, save_as_mrc, read_gemmi_map, compare_gemmi_grids
    from emmer.ndimage.map_tools import get_center_of_mass
    from emmer.ndimage.map_tools import compute_real_space_correlation
    
    
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
                                          return_grid=True, align_output=True, verbose=True, set_refmac_blur=True, blur=add_blur)
    
    ## Normalise intensity levels of refmap_data
    refmap_data_normalised = normalise_intensity_levels(
        from_emmap=refmap_data, to_levels=[emmap_data.min(),emmap_data.max()])
    
    ## Output filename
    reference_map_path = model_path[:-4]+"_4locscale.mrc"
    save_as_mrc(map_data=refmap_data_normalised,output_filename=reference_map_path, apix=grid_simulated.spacing, origin=0)   
    
    ## Checklist: 
    
    correlation = compute_real_space_correlation(emmap_path, reference_map_path)
    
    grid_comparison = compare_gemmi_grids(read_gemmi_map(emmap_path, return_grid=True)[1], 
                                          read_gemmi_map(reference_map_path,return_grid=True)[1])
    # Since grid comparison and correlation are critical, they are done on the saved filesystem 
    # and not on the files in memory. This is to avoid any errors that might have happened during 
    # save_as_mrc operation
    
    center_of_mass_experimental = get_center_of_mass(emmap_data,apix=grid_input.spacing)
    center_of_mass_simulated = get_center_of_mass(refmap_data_normalised,apix=grid_simulated.spacing)
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
    
def run_mapmask(emmap_path, return_same_path=False):
    '''
    Function to generate a XYZ corrected output using CCPEM-Mapmask tool

    Parameters
    ----------
    emmap_path : str
        

    Returns
    -------
    mapmasked_path : str

    '''
    import os
    from subprocess import run
    path_to_locscale = check_dependencies()['locscale']
    
    mapmask_bash_script = path_to_locscale + "/locscale/utils/mapmask.sh"
    
    if return_same_path:
        xyz_output_map = emmap_path
    else:
        xyz_output_map = "/".join(emmap_path.split('/')[:-1]+   ["xyz_"+emmap_path.split(sep='/')[-1]])
    run([mapmask_bash_script,emmap_path,xyz_output_map])
    
    return xyz_output_map
    
