#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 17 11:31:51 2020

@author: alok
"""

import numpy as np
from locscale.utils.plot_tools import tab_print

tabbed_print = tab_print(2)
tprint = tabbed_print.tprint
def prepare_sharpen_map(emmap_path,wilson_cutoff,fsc_resolution,add_blur=20,return_processed_files=False, output_file_path=None,verbose=True):
    '''
    Function to prepare target map for refinement 

    Parameters
    ----------
    emmap_path : str
        Path to the EM map which needs sharpening. Example: 'path/to/map.mrc'
    wilson_cutoff : float
        Wilson cutoff for determining the B factor of the EM map
    fsc_resolution : float
        Resolution at which FSC curve is computed used to compute B factor
    add_blur : float, optional
        Additional blur to be added to the EM map. The default is 20   
    
    Returns
    -------
    sharpen_map_path : str
        Path to the sharpened EM map. Example: 'path/to/sharpen_map.mrc'


    '''

    from locscale.include.emmer.ndimage.profile_tools import compute_radial_profile, estimate_bfactor_through_pwlf, frequency_array
    from locscale.include.emmer.ndimage.map_utils import average_voxel_size, save_as_mrc
    from locscale.include.emmer.ndimage.map_tools import sharpen_maps, estimate_global_bfactor_map
    from locscale.include.emmer.ndimage.filter import apply_filter_to_map
    import mrcfile
    import warnings

    emmap_mrc = mrcfile.open(emmap_path)
    emmap_unsharpened = emmap_mrc.data
    apix=average_voxel_size(emmap_mrc.voxel_size)
    
    if verbose:
        tprint("Estimating global B-factor using breakpoints\n")
    
    ### Estimate global B-factor ### 
    
    bfactor, z, slopes, fit = estimate_global_bfactor_map(emmap=emmap_unsharpened, apix=apix,
                                    wilson_cutoff=wilson_cutoff, fsc_cutoff=fsc_resolution)
    
    if verbose:
        tprint("Global bfactor: {:.3f}\n".format(-1*bfactor))
        tprint("Breakpoints: {}".format(1/np.sqrt(z)))
        tprint("Slopes: {}".format(slopes))
    
    if verbose:
        tprint("Calculating final b-factor of target map")
    if add_blur != 0:
        bfactor  += add_blur  ## Use add_blur if you wanna add blur to the emmap before refining
    if verbose:
        tprint("Final overall bfactor of emmap expected to be {:.2f}".format(add_blur))
    

    ### Globally sharpen EM Map ###
    sharpened_map = sharpen_maps(emmap_unsharpened, apix=apix, global_bfactor=bfactor)
       
    bfactor_final, z_final, slopes_final, fit_final = estimate_global_bfactor_map(emmap=sharpened_map, apix=apix, 
                                                    wilson_cutoff=wilson_cutoff, fsc_cutoff=fsc_resolution)
    
    if verbose:
        tprint("Final overall bfactor of emmap computed to be {:.2f}".format(-1*bfactor_final))
        tprint("Breakpoints: {}".format(1/np.sqrt(z_final)))
        tprint("Slopes: {}".format(slopes_final))
        
    
    output_filename = emmap_path[:-4] +"_global_sharpened.mrc"
    output_filename_filtered_map = emmap_path[:-4] +"_global_sharpened_filtered.mrc"
    
    save_as_mrc(map_data=sharpened_map, output_filename=output_filename, apix=apix)
    apply_filter_to_map(output_filename, dmin=fsc_resolution, output_filename=output_filename_filtered_map)
    
    if return_processed_files:
        if verbose:
            tprint("Returning: sharpend_map_path (filtered at FSC), [rp_unsharp, rp_sharp, bfactor]")
        return output_filename_filtered_map, fit
    else:                    
        return output_filename_filtered_map

def run_FDR(emmap_path,window_size,fdr=0.01,verbose=True,filter_cutoff=None):
    '''
    Function to calculate FDR mask for a map

    Parameters
    ----------
    emmap_path : string
        Path to the EM Map which needs thresholding. Example: 'path/to/map.mrc'
    window_size : int
        Window size required for FDR thresholding
    fdr : float, optional
        FDR value to be used for thresholding. The default is 0.01
    verbose : bool, optional
        Print statistics if True. The default is True.

    Returns
    -------
    mask_path : string
        path to mask file. Example: 'path/to/mask.mrc'

    '''
    import os, sys
    from subprocess import run, PIPE
    import mrcfile
    # Preprocessing EM Map Path
    
    # Apply filter if filter_cutoff is not None

    if verbose:
        tprint("Now starting FDR procedure using the following parameters: \n"
                 "Window size: "+str(window_size)+"\n"
                 "Filter cutoff: "+str(filter_cutoff))
    
        
    from locscale.include.emmer.ndimage.map_utils import average_voxel_size, compute_FDR_confidenceMap_easy, save_as_mrc
    emmap = mrcfile.open(emmap_path).data
    voxel_size_record = mrcfile.open(emmap_path).voxel_size
    apix = average_voxel_size(voxel_size_record)
    fdr = fdr
    fdr_mask, fdr_threshold = compute_FDR_confidenceMap_easy(
        emmap, apix=apix, fdr=fdr, window_size=window_size, 
        lowPassFilter_resolution=filter_cutoff, remove_temp_files=False)
    
    
    emmap_path_without_ext = emmap_path[:-4]
    mask_path = emmap_path_without_ext + "_confidenceMap.mrc"
    
    save_as_mrc(fdr_mask, output_filename=mask_path, 
                apix=voxel_size_record.tolist(), origin=0)
    
    if os.path.exists(mask_path):
        if verbose:
            tprint("FDR Procedure completed. \n"+
                    "Mask path: "+mask_path+"\n")
        
        return mask_path
    else:
        tprint("FDR process failed. Returning none")
        return None
         
def run_pam(emmap_path,mask_path,threshold,num_atoms,method,bl,
            g=None,friction=None,scale_map=None,scale_lj=None,total_iterations=100,verbose=True):
    '''
    Function to build a pseudo atomic model     

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
    from locscale.preprocessing.pseudomodel_solvers import main_solver3D, main_solver_kick
    from locscale.preprocessing.pseudomodel_classes import extract_model_from_mask
    
    #### Input ####

    mrc = mrcfile.open(emmap_path)
    emmap = mrc.data
    voxelsize = mrc.voxel_size.x

    mask = mrcfile.open(mask_path).data

    ### Extract a random pseudo-atomic model from the masked map ###
    pseudomodel = extract_model_from_mask(mask,num_atoms,threshold=threshold)
    
    emmap_shape = emmap.shape
    unitcell = gemmi.UnitCell(emmap_shape[0]*voxelsize,emmap_shape[1]*voxelsize,emmap_shape[2]*voxelsize,90,90,90)
    
    outputlogfilepath = os.path.join(os.path.dirname(emmap_path),"pseudomodel_log.txt")
    output_file = open(outputlogfilepath,"w")

    ### Run the solver depending on user input: gradient or random ###
    if verbose:
        tprint("Running pseudoatomic model generator to add "+str(num_atoms)+" atoms inside the volume using the method: "+method)
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
            emmap,gx,gy,gz,pseudomodel,g=g,friction=friction,min_dist_in_angst=bl,voxelsize=voxelsize,dt=0.1,
            capmagnitude_lj=100,epsilon=1,scale_lj=scale_lj,capmagnitude_map=100,scale_map=scale_map,
            total_iterations=total_iterations, compute_map=None,emmap_path=None,mask_path=None,
            returnPointsOnly=True,integration='verlet',verbose=False, myoutput=output_file)

        mask_name = mask_path[:-4]
        pseudomodel_path = mask_name+"_gradient_pseudomodel.pdb"

    elif method=='random' or method=='kick' or method == 'random_placement_with_kick':
        arranged_points = main_solver_kick(
                pseudomodel,min_dist_in_angst=bl,voxelsize=voxelsize,total_iterations=total_iterations,
                returnPointsOnly=True,verbose=False, myoutput=output_file)
        mask_name = mask_path[:-4]
        pseudomodel_path = mask_name+"_kick_pseudomodel.pdb"
    

    ### Save the pseudo-atomic model ###    
    arranged_points.write_pdb(pseudomodel_path,voxelsize=voxelsize,unitcell=unitcell)
    
    
    if os.path.exists(pseudomodel_path):    
        if verbose:
            tprint("The location of the pseudomodel generated is: "+pseudomodel_path+'\n\n')
        return pseudomodel_path
    else:
        tprint("uhhu, something wrong with the pseudomodel generator! Returning None")        
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
        return True
    else:
        return False
    
        
    
def run_refmac_servalcat(model_path, map_path,resolution,  num_iter, only_bfactor_refinement, refmac5_path=None, verbose=True):
    '''
    Function to run Refmac to refine the model and generate a new model with refined B-factors.

    Parameters
    ----------
    model_path : string
        Path of the model to be refined
    map_path : string
        Path of the map to be used for refinement
    resolution : float
        Resolution of the map
    num_iter : int
        Number of iterations to run Refmac
    only_bfactor_refinement : bool
        If True, only the B-factors will be refined. If False, the model will be refined using the map.
    refmac5_path : string
        Path of the refmac5 executable. If None, the default path will be used.
    verbose : bool
        If True, the output of Refmac will be printed.
    
    Returns
    -------
    refined_model_path : string
        Path of the refined model
    
    '''
    
    import os
    from subprocess import run, PIPE, Popen
    from locscale.include.emmer.pdb.pdb_utils import get_bfactors, set_atomic_bfactors
    import mrcfile
    
    #### Input ####
    model_name = os.path.basename(model_path)
    servalcat_uniform_bfactor_input_path = model_path[:-4]+"_uniform_biso.pdb"

    #### Set the starting bfactor of the atomic model to 40 before refinement ####

    set_atomic_bfactors(in_model_path=model_path, b_iso=40, out_file_path=servalcat_uniform_bfactor_input_path)
    
    ### Run Servalcat after preparing inputs ###
    output_prefix = model_name[:-4]+"_servalcat_refined"
    servalcat_command = ["servalcat","refine_spa","--model",servalcat_uniform_bfactor_input_path,"--resolution",str(round(resolution, 2)), "--map", map_path, "--ncycle",str(int(num_iter)), "--output_prefix",output_prefix]
    if only_bfactor_refinement:
        servalcat_command += ["--keywords","refi bonly","refi type unre"]

    if refmac5_path is not None:
        servalcat_command += ["--exe",refmac5_path]

    if verbose:       
        tprint("Command line for servalcat: \n")    
        tprint(" ".join(servalcat_command))
    
    ## Create log file for servalcat
    servalcat_log_path = os.path.join(os.path.dirname(model_path),"servalcat_log.txt")
    servalcat_log_file = open(servalcat_log_path,"w")
    
    refmac_output = run(servalcat_command, stdout=servalcat_log_file)
    refined_model_path = output_prefix+".pdb"
    bfactors = get_bfactors(in_model_path=refined_model_path)
    
    if os.path.exists(refined_model_path):
        if verbose: 
            tprint("The refined PDB model is: "+refined_model_path+"\n\n")    
            tprint("B factor range: \t ({:.2f} to {:.2f})".format(min(bfactors),max(bfactors)))
            
        return refined_model_path
    else:
        tprint("Uhhoh, something wrong with the REFMAC procedure. Returning None")
        return None


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
    from locscale.include.emmer.pdb.pdb_to_map import pdb2map
    from locscale.include.emmer.ndimage.map_utils import average_voxel_size, save_as_mrc, read_gemmi_map, compare_gemmi_grids
    from locscale.include.emmer.ndimage.map_tools import get_center_of_mass
    from locscale.include.emmer.ndimage.map_tools import compute_real_space_correlation
    
    
    if verbose: 
        tprint("Now simulating Reference Map using Refined Atomic Model")
    
    # Read inputs from filesystem
    
    emmap_data, grid_input = read_gemmi_map(emmap_path, return_grid=True)
    mask = read_gemmi_map(mask_path)
    pdb_structure = gemmi.read_structure(model_path)
    
    
    ## Generate parameters of the simulated map
    voxelsize = grid_input.spacing   
    unitcell = grid_input.unit_cell
    
    ## Simulate a reference map from the input atomic model in the pdb_structure variable
    
    refmap_data, grid_simulated = pdb2map(input_pdb=pdb_structure, unitcell=unitcell, size=emmap_data.shape,
                                          return_grid=True, align_output=True, verbose=False, set_refmac_blur=True, blur=add_blur)
    

    refmap_data_normalised = refmap_data
    
    
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
            tprint("The reference map is at: "+reference_map_path+"\n\n")
            for key, value in reporter.items():
                tprint(key+":\t"+str(value))
        return reference_map_path
    else:
        tprint("Reference map was not generated. Returning none")
        return None
    

def check_axis_order(emmap_path, return_same_path=False):
    '''
    Function to generate a XYZ corrected output using Gemmi

    Parameters
    ----------
    emmap_path : str
        

    Returns
    -------
    mapmasked_path : str

    '''
    import os
    from subprocess import run, PIPE
    from locscale.include.emmer.ndimage.map_utils import read_gemmi_map, ZYX_to_XYZ, save_as_mrc
           
    emmap, grid = read_gemmi_map(emmap_path, return_grid=True)

    ## Check if the map is in the right order
    if grid.axis_order.name == "XYZ":
        xyz_emmap_path = emmap_path
    elif grid.axis_order.name == "ZYX":
        ## Flip and rotate the map
        xyz_emmap = ZYX_to_XYZ(emmap)
        xyz_emmap_path = os.path.join(os.path.dirname(emmap_path), "xyz_"+os.path.basename(emmap_path))
        save_as_mrc(map_data=xyz_emmap,output_filename=xyz_emmap_path, apix=grid.spacing)
    else:
        print("### Warning: Map is not in the right order. Please check the axis order of the map")
        print("Axis order of the map: "+grid.axis_order.name)
        print("Using the same map as input")
        xyz_emmap_path = emmap_path
    
    return xyz_emmap_path
    


############# CODE HELL #############

# def run_mapmask(emmap_path, return_same_path=False):
#     '''
#     Function to generate a XYZ corrected output using CCP4-Mapmask tool
#     DEPRECATED: Use check_axis_order instead

#     Parameters
#     ----------
#     emmap_path : str
        

#     Returns
#     -------
#     mapmasked_path : str

#     '''
#     import os
#     from subprocess import run, PIPE
#     from locscale.utils.file_tools import get_locscale_path
#     path_to_locscale = get_locscale_path()
    
#     mapmask_bash_script = os.path.join(path_to_locscale , "locscale","utils","mapmask.sh")
    
#     if return_same_path:
#         xyz_output_map = emmap_path
#     else:
#         xyz_output_map = "/".join(emmap_path.split('/')[:-1]+   ["xyz_"+emmap_path.split(sep='/')[-1]])
#     run([mapmask_bash_script,emmap_path,xyz_output_map], stdout=PIPE)
    
#     return xyz_output_map


# """ def run_refmac(model_path,map_path,resolution,  num_iter,only_bfactor_refinement,verbose=True):
#     import os
#     from subprocess import run, PIPE
#     from locscale.include.emmer.pdb.pdb_utils import get_bfactors
#     import mrcfile
    
#     path_to_locscale = get_locscale_path()
#     path_to_ccpem = get_locscale_path()['ccpem']
#     path_to_ccp4 = get_locscale_path()['ccp4']
    
#     if only_bfactor_refinement:
#         path_to_run_refmac = os.path.join(path_to_locscale,"locscale","utils","run_refmac.sh")
#     else:
#         path_to_run_refmac = os.path.join(path_to_locscale,"locscale","utils","run_refmac_restrained.sh")
        
    
#     model_name = model_path[:-4]
#     emmap_mrc = mrcfile.open(map_path)
#     map_dims = emmap_mrc.header.cella.tolist()
    
#     refmac_command_line = "bash "+path_to_run_refmac+" "+model_path+" "+model_name+" "+map_path+" "+str(round(resolution,2))+" "+path_to_ccpem+" "+path_to_ccp4+" "+str(map_dims[0])+" "+str(map_dims[1])+" "+str(map_dims[2])+" "+str(num_iter)
    
    
#     if verbose:
#         print("Running REFMAC to refine the pseudo-atomic model using \n"+
#               "Path to run_refmac: "+path_to_run_refmac+"\n"+
#               "Command line: \n"+refmac_command_line)
        
#     refmac_output = run(refmac_command_line.split())
#     refined_model_path = model_name+"_refmac_refined.pdb"
#     bfactors = get_bfactors(in_model_path=refined_model_path)
    
#     if os.path.exists(refined_model_path):
#         if verbose: 
#             print("The refined PDB model is: "+refined_model_path+"\n\n")    
#             print("B factor range: \t ({:.2f} to {:.2f})".format(min(bfactors),max(bfactors)))
            
#         return refined_model_path
#     else:
#         print("Uhhoh, something wrong with the REFMAC procedure. Returning None")
#         return None
#      """