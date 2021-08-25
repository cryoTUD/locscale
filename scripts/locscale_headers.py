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
from pdb_tools import *

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


def run_FDR(emmap_path,window_size,verbose=True,filter_cutoff=None):
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
    
    # Preprocessing EM Map Path
    
    # Apply filter if filter_cutoff is not None
    if filter_cutoff is not None:
        filtered_emmap = apply_filter_to_map(emmap_path, filter_cutoff)
        emmap_path = filtered_emmap
        
    path_to_FDR_script = path_to_ccpem+"/lib/py2/FDRcontrol.pyc"
    fdr_command_line = "ccpem-python "+path_to_FDR_script+" --em_map "+emmap_path+" -method BY --testProc rightSided --window_size "+str(window_size)
    
    statistics = {}
    statistics['window_size'] = window_size
    statistics['fdr_script'] = path_to_FDR_script
    statistics['fdr_command_line'] = fdr_command_line
    
    if verbose:
        print("Now starting FDR procedure using the following parameters: \n"
              "Window size: "+str(window_size)+"\n"
              "FDR Script Location: \n"+path_to_FDR_script+
              "FDR Command Line: \n"+fdr_command_line)
        
        
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
        

def measure_mask_parameters(mask_path,edge_threshold=1,protein_density=1.35,average_atomic_weight=13.14,verbose=True,detailed_report=False):
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
        Found using 54% carbon, 20% oxygen and 16% nitrogen. The default is 12.066.
    verbose : bool, optional
        Print statistics if True. The default is True.

    Returns
    -------
    num_atoms : int
        Estimated number of atoms based on mask volume, protein density and average atomic weight
    

    '''
    mask_mrc = mrcfile.open(mask_path)
    voxelsize = mask_mrc.voxel_size.x
    ang_to_cm = 1e-8
    mask = mask_mrc.data
    mask.setflags(write=1)
    mask[mask<edge_threshold] = 0
    mask[mask>=edge_threshold] = 1
    mask_vol = mask.sum()*(voxelsize*ang_to_cm)**3
    mask_vol_A3 = mask.sum()*voxelsize**3
    #print("\n Volume of the mask generated is: "+str(mask.sum())+" A$^3$ \n")
    # Calculate number of atoms
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

def run_pam_using_phenix(emmap_path,mask_path,threshold,num_atoms,method,g,bl,friction,scale_map,scale_lj,
                         total_iterations,use_phenix,resolution,myoutput):
    print("************************************************")
    print(datetime.now())
    if use_phenix is False:
        mrc = mrcfile.open(emmap_path)
        emmap = mrc.data
        emmap = emmap.clip(min=0)
        voxelsize = mrc.voxel_size.x
        gz,gy,gx = np.gradient(emmap)
        
        mask = mrcfile.open(mask_path).data
        pseudomodel = extract_model_from_mask(mask,num_atoms,threshold=threshold)
        
        print("Next, I'll add "+str(num_atoms)+" empty atoms at random points throughout the mask to generate a pseudo-atomic model. \n")
        print("************************************************")
        if method=='gradient':
            arranged_points = main_solver3D(
                emmap,gx,gy,gz,pseudomodel,g=g,friction=friction,min_dist_in_angst=bl,voxelsize=voxelsize,dt=0.1,capmagnitude_lj=100,epsilon=1,scale_lj=scale_lj,
                capmagnitude_map=100,scale_map=scale_map,total_iterations=total_iterations, path_for_gemmi_models=None,emmap_path=None,mask_path=None,returnPointsOnly=True,
                verbose=True,integration='verlet',myoutput=myoutput
                )
            mask_name = mask_path[:-4]
            pseudomodel_path = mask_name+"_gradient_pseudomodel.pdb"
            arranged_points.write_pdb(pseudomodel_path)
            #emmi_model_gradient = convert_to_gemmi_model(arranged_points.list,voxelsize=voxelsize)
            #rite_pdb(gemmi_model_gradient,pseudomodel_path) 
        
        elif method=='kick':
            arranged_points = main_solver_kick(
                    pseudomodel,min_dist_in_angst=1.8,voxelsize=voxelsize,total_iterations=99,returnPointsOnly=True,verbose=True)
            mask_name = mask_path[:-4]
            pseudomodel_path = mask_name+"_kick_pseudomodel.pdb"
            arranged_points.write_pdb(pseudomodel_path)
            
            #gemmi_model_kick = convert_to_gemmi_model(arranged_points.list,voxelsize=voxelsize)
            #write_pdb(gemmi_model_kick,pseudomodel_path)
    else:
        # Generate helix_strands pdb file
        # Covnert map to MTZ
        
        mask = mrcfile.open(mask_path).data
        voxelsize = mrcfile.open(mask_path).voxel_size.x
        mask_name = mask_path[:-4]
        
        mtz_filename = emmap_path[:-4]+".mtz"
        command_line_map_to_mtz = "phenix.map_to_structure_factors "+emmap_path+" d_min="+str(resolution)+" output_file_name="+mtz_filename
        print("Now running: \n")
        print(command_line_map_to_mtz+'\n')
        myoutput.write('\nPhenix Command Line: \n')
        myoutput.write(command_line_map_to_mtz+'\n')
        mtzoutput = run(command_line_map_to_mtz.split(),stdout=myoutput)
        
        # Convert MTZ to PDB
        hs_pdb_filename = emmap_path[:-4]+"_helix_strand.pdb"
        command_line_hs = "phenix.find_helices_strands "+mtz_filename+" output_model="+hs_pdb_filename
        print("Now running: \n")
        print(command_line_hs+'\n')
        myoutput.write('\nPhenix Command Line: \n')
        myoutput.write(command_line_hs+'\n')
        mtzoutput = run(command_line_hs.split(),stdout=myoutput)
        
        print("Helices and Strands found. Now adding remaining atoms at random locations")
        # Now add rest of atoms
        helix_strand_model = get_model_from_gemmi_pdb(hs_pdb_filename,emmap_path)
        num_helix_strand_atoms = len(helix_strand_model.list)
        num_remaining_atoms = num_atoms - num_helix_strand_atoms
        
        remaining_model = extract_model_from_mask(mask,num_remaining_atoms,threshold=1,ignore_these=helix_strand_model.extract_mrc_positions())
        for atom in remaining_model.list:
             atom.pdb_position = atom.position.scale(voxelsize)
        pseudomodel_path = mask_name+"_pseudomodel.pdb"
        
        remaining_model.combine(helix_strand_model)
        remaining_model.write_pdb(pseudomodel_path)
        #gemmi_model = convert_to_gemmi_model(remaining_model.list,voxelsize=voxelsize)
        #write_pdb(gemmi_model,pseudomodel_path) 

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
        
        gz,gy,gx = np.gradient(emmap)
        arranged_points = main_solver3D(
            emmap,gx,gy,gz,pseudomodel,g=g,friction=friction,min_dist_in_angst=bl,voxelsize=voxelsize,dt=0.1,capmagnitude_lj=100,epsilon=1,scale_lj=scale_lj,
            capmagnitude_map=100,scale_map=scale_map,total_iterations=total_iterations, path_for_gemmi_models=None,emmap_path=None,mask_path=None,returnPointsOnly=True,
            integration='verlet',verbose=True)
        mask_name = mask_path[:-4]
        pseudomodel_path = mask_name+"_gradient_pseudomodel.pdb"

    elif method=='random_placement_with_kick':
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
    
    
    
def run_refmap(model_path,emmap_path,mask_path,verbose=True):
    if verbose: 
        print("Now simulating Reference Map using Refined Atomic Model")
    reference_map = model_path[:-4]+"_4locscale.mrc"
    emmap_mrc = mrcfile.open(emmap_path)
    voxelsize = emmap_mrc.voxel_size.x
    
    mask_data = mrcfile.open(mask_path).data
    refmap_data = pdb_to_map(pdb_path=model_path,vsize=voxelsize,set_zero_origin_and_crop=True,save_mrc_path=reference_map,remove_waters=False)
    emmap_data = emmap_mrc.data
    
    
    if os.path.exists(reference_map):
        if verbose: 
            correlation = compute_real_space_correlation(emmap_data, refmap_data)
            masked_correlation = compute_real_space_correlation(emmap_data*mask_data, refmap_data*mask_data)
            print("The reference map is at: "+reference_map+"\n\n")
            print("MAP STATISTICS: \n"+
                  "Real Space Correlation with EM Map: "+str(round(correlation,3))+"\n"+
                  "Real Space Correlation with EM Map (masked): "+str(round(masked_correlation,3))+"\n")
            
        return reference_map
    else:
        print("Uhhu, something wrong with the Reference Map generation. Returning none")
        return None
    
def run_refmap2(model_path,emmap_path,mask_path,verbose):
    print("************************************************")
    print(datetime.now())
    print("Just one more step.. I shall now use the refined b factors from the Refmac step to generate a Reference Map, which can be used for locscale at this point!")
    print("************************************************")
    path_to_refmap_script = path_to_locscale+"/source/prepare_locscale_input.py"
    refmap_command_line = "phenix.python "+path_to_refmap_script+" -mc "+model_path+" -em "+emmap_path+" -ma "+mask_path
    print("Now running: \n")
    print(refmap_command_line+'\n')

    refmap_output = run(refmap_command_line.split(),stdout=PIPE)
    reference_map = model_path[:-4]+"_4locscale.mrc"
    new_emmap = emmap_path[:-4]+"_4locscale.mrc"
    new_mask_path = mask_path[:-4]+"_4locscale.mrc"
    
    if os.path.exists(reference_map):
        print("The reference map is at: "+reference_map+"\n\n")
        return reference_map, new_emmap, new_mask_path
    else:
        print("Uhhu, something wrong with the Reference Map generation. Please check the log file")
        return None

def run_mapmask(map_path):
    mapmask_file =path_to_locscale +"/scripts/mapmask.sh"
    
    command_line = "bash "+mapmask_file+" "+map_path
    print(command_line)
    output =run(command_line.split(),stdout=PIPE)
    output_map = "xyz_"+map_path
    if os.path.exists(output_map):
        print("MAPMASK succesful")
        return output_map
    else:
        print("MAPMASK unsuccessful")
        return None
def run_locscale(emmap_path,refmap_path,apix,wsize,myoutput):
    print("************************************************")
    print(datetime.now())
    print("Now for the final step: Locscale! \n ")
    print("************************************************")
    path_to_run_locscale = path_to_locscale+"/scripts/run_locscale.zsh"
    directory = '/'.join(emmap_path.split(sep='/')[:-1])
    emmap_name = emmap_path.split(sep='/')[-1]
    refmap_name = refmap_path.split(sep='/')[-1]
    xyz_emmap_path = directory+'xyz_'+emmap_name
    xyz_refmap_path = directory+'xyz_'+refmap_name
    locscale_command_line = "zsh "+path_to_run_locscale+" "+emmap_path+" "+xyz_emmap_path+" "+refmap_path+" "+xyz_refmap_path+" "+str(apix)+" "+str(wsize)+" "+path_to_ccpem+" "+path_to_ccp4
    print("Now running: \n")
    print(locscale_command_line+'\n')    
    myoutput.write('\nLocscale Command Line: \n')
    myoutput.write(locscale_command_line+'\n')
    run(locscale_command_line.split(),stdout=myoutput)


def launch_locscalev02(args):
    print("***** Welkom to Locscale version 0.2! ***** \n Now you can sharpen your boring and blurry EM density maps *WITHOUT* needing a reference map to start with! I hope you selected an unfiltered EM map while runnning me!\n\n")
    emmap_path = args.em_map
    method = args.method
    print("You can find the detailed verbose outout in the file locscale_output.txt stored in this same folder. ")
    directory = '/'.join(emmap_path.split(sep='/')[:-1])
    emmap_name = emmap_path[:-4]
    output_filepath = directory+emmap_name+'_locscale_output.txt'
    myoutput = open(output_filepath,'w+')
    #sys.stdout = myoutput

    resolution = round(float(args.resolution),2)
    wsize_fdr = int(args.wsize_fdr)
    wsize_locscale = int(args.wsize_locscale)
    g = int(args.gradient_scale)
    bl = float(args.bond_length)
    friction = float(args.friction)
    scale_map = float(args.scale_map)
    scale_lj = float(args.scale_lj)
    only_pseudo = bool(args.only_pseudo)
    total_iterations = int(args.total_iterations)
    mask_threshold = float(args.threshold)
    use_phenix = bool(args.use_phenix)
    
    if not only_pseudo:
        mask_path,mask_mrc = run_FDR(emmap_path,wsize_fdr,myoutput)
        num_atoms,voxelsize,maskdims = measure_mask_parameters(args,mask_mrc,mask_threshold)
        pseudomodel_path = run_pam(emmap_path,mask_path,mask_threshold,num_atoms,method,g,bl,friction,scale_map,scale_lj,total_iterations,use_phenix,resolution,myoutput)
    
    if only_pseudo:
        if args.mask_path is not None:
            mask_path = args.mask_path
            mask_mrc = mrcfile.open(mask_path)
            num_atoms,voxelsize,maskdims = measure_mask_parameters(args,mask_mrc,mask_threshold)
            pseudomodel_path = run_pam(emmap_path,mask_path,mask_threshold,num_atoms,method,g,bl,friction,scale_map,scale_lj,total_iterations,use_phenix,resolution,myoutput)
            
        elif args.num_atoms is not None:
            num_atoms = int(args.nrandom_placeum_atoms)
            mask_path = create_pseudo_mask(emmap_path,float(args.threshold))
            voxelsize = mrcfile.open(emmap_path).voxel_size.x
            mask = mrcfile.open(mask_path).data
            maskshape = mask.shape
            maskdims = [maskshape[0]*voxelsize,maskshape[1]*voxelsize,maskshape[2]*voxelsize]
            pseudomodel_path = run_pam(emmap_path,mask_path,0.99,num_atoms,method,g,bl,friction,scale_map,scale_lj,total_iterations,use_phenix,resolution,myoutput)
    
    
    
    if not only_pseudo:
        refined_model_path = run_refmac(pseudomodel_path,pseudomodel_path[:-4],emmap_path,resolution,maskdims,myoutput)
        
        reference_map_path = run_refmap(refined_model_path,emmap_path,mask_path,myoutput)
        emmap_path_4locscale = emmap_path[:-4]+'_4locscale.mrc'
        
        run_locscale(emmap_path_4locscale,reference_map_path,voxelsize,wsize_locscale,myoutput)
    
    if os.path.exists('loc_scale.mrc'):
        print('Locscale procedure completed. You can find your scaled map here >>>> '+os.getcwd()+'/loc_scale.mrc')
    
    print("********************* THE END *****************************")
    myoutput.close()
