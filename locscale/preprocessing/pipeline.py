
def get_modmap(modmap_args):
    '''
    Function to generate a model map using pseudo-atomic model

    Parameters
    ----------
    modmap_args : dict
    Returns
    -------
    pseudomodel_modmap : str
        path/to/modmap.mrc

    '''
    from locscale.preprocessing.headers import run_FDR, run_pam, run_refmac_servalcat, run_refmap, prepare_sharpen_map, is_pseudomodel, run_servalcat_iterative
    from locscale.include.emmer.ndimage.map_utils import measure_mask_parameters, average_voxel_size
    from locscale.include.emmer.ndimage.map_tools import estimate_global_bfactor_map, estimate_global_bfactor_map_standard
    from locscale.include.emmer.pdb.pdb_tools import find_wilson_cutoff, add_pseudoatoms_to_input_pdb
    from locscale.include.emmer.pdb.pdb_utils import get_bfactors, add_atomic_bfactors, shift_bfactors_by_probability
    from locscale.utils.plot_tools import tab_print
    import mrcfile
    from scipy.stats import norm
    import pickle
    import numpy as np
    import os
    
    ###########################################################################
    # Extract the inputs from the dictionary
    ###########################################################################
    tabbed_print = tab_print(2)
    
    emmap_path = modmap_args['emmap_path']
    halfmap_paths = modmap_args['halfmap_paths']
    mask_path = modmap_args['mask_path_raw']
    pdb_path = modmap_args['pdb_path']
    pseudomodel_method = modmap_args['pseudomodel_method']
    pam_distance = modmap_args['pam_distance']
    pam_iteration = modmap_args['pam_iteration']
    fsc_resolution = modmap_args['fsc_resolution']
    refmac_iter = modmap_args['refmac_iter']
    add_blur = modmap_args['add_blur']
    skip_refine = modmap_args['skip_refine']
    refmac5_path = modmap_args['refmac5_path']
    pg_symmetry = modmap_args['pg_symmetry']
    model_resolution = modmap_args['model_resolution']
    molecular_weight = modmap_args['molecular_weight']
    build_ca_only = modmap_args['build_ca_only']
    verbose = modmap_args['verbose']
    #Cref = modmap_args['Cref']
    complete_model = modmap_args['complete_model']
    averaging_window = modmap_args['averaging_window']
    mask_threshold = modmap_args['mask_threshold']
    cif_info = modmap_args['cif_info']

    if verbose:
        print("."*80)
        print("Running model-map generation pipeline \n")


    if verbose:
        tabbed_print.tprint("Model map arguments: \n")
        ## Print keys and values of dictionary in a nice format
        for key, value in modmap_args.items():
            if key == "Cref":
                # Print Cref shape
                if value is not None:
                    tabbed_print.tprint("{} : {}".format(key, value.shape))
                else:
                    tabbed_print.tprint("{} : {}".format(key, value))
            else:
                tabbed_print.tprint("{:<20} : {}".format(key, value))

    #########################################################################
    # Open data files and collect required inputs
    # #######################################################################    
    emmap_mrc = mrcfile.open(emmap_path)
    apix = average_voxel_size(emmap_mrc.voxel_size)
    
    pam_bond_length = pam_distance
    pam_method = pseudomodel_method
    pam_iteration = pam_iteration
    resolution = fsc_resolution
    verbose = verbose
    ###########################################################################
    # Stage 1: Check the required number of atoms for the pseudomodel
    ###########################################################################
    if molecular_weight is None:
        num_atoms,mask_dims = measure_mask_parameters(mask_path,verbose=False, edge_threshold=mask_threshold)
    else:
        avg_mass_per_atom = 13.14  #amu
        num_atoms = int(molecular_weight * 1000.0 / avg_mass_per_atom)
    ###########################################################################
    # Stage 1a: Check if the user requires to build only Ca atoms
    ###########################################################################
    if build_ca_only:
        num_atoms = int(num_atoms/9)  ## Assuming 9 atoms per residue
        pam_bond_length = 3.8  ## Ca atom distances for secondary structures
        pam_method = 'gradient'  ## use this exclusively for Gradient
        if pam_method != 'gradient':
            tabbed_print.tprint("Using gradient method for building pseudo-atomic model!\
                Not using user input:\t {}".format(pam_method))
    ###########################################################################
    # Stage 1b: If user has not provided a PDB path then build a 
    # pseudomodel using the run_pam() routine else use the PDB path directly
    ###########################################################################
    if pdb_path is None:
        if verbose:
            print("."*80)
            print("You have not entered a PDB path, running pseudo-atomic model generator!")
        input_pdb_path = run_pam(emmap_path=emmap_path, mask_path=mask_path, threshold=mask_threshold, num_atoms=num_atoms, 
                                   method=pam_method, bl=pam_bond_length,total_iterations=pam_iteration,verbose=verbose)
        if input_pdb_path is None:
            print("Problem running pseudo-atomic model generator. Returning None")
            return None
        final_chain_counts = None
    else:
        
        if complete_model:
            if verbose:
                print("."*80)
                print("Adding pseudo-atoms to the regions of the mask that are not modelled by the user-provided PDB")
            integrated_structure, final_chain_counts, difference_mask_path = add_pseudoatoms_to_input_pdb(
                pdb_path=pdb_path, mask_path=mask_path, emmap_path=emmap_path,\
                averaging_window=averaging_window, pseudomodel_method=pam_method, pseudomodel_iteration=pam_iteration, \
                mask_threshold=mask_threshold, fsc_resolution=fsc_resolution, \
                return_chain_counts=True, return_difference_mask=True) 

            input_pdb_path = pdb_path[:-4] + '_integrated_pseudoatoms.pdb'
            integrated_structure.write_pdb(input_pdb_path)
        else:
            final_chain_counts = None
            if verbose:
                print("."*80)
                print("Using user-provided PDB path: {}".format(pdb_path))    
            input_pdb_path = pdb_path
    ###########################################################################
    # Stage 2: Refine the reference model usign servalcat
    ###########################################################################
    if is_pseudomodel(input_pdb_path):
        pseudomodel_refinement = True
    else:
        pseudomodel_refinement = False
            
    wilson_cutoff = find_wilson_cutoff(mask_path=mask_path, return_as_frequency=False, verbose=False)
    
    # #############################################################################
    # # Stage 2a: Prepare the target map for refinement by globally sharpening
    # # the input map
    # #############################################################################
    # if verbose:
    #     print("."*80)
    #     print("Preparing target map for refinement\n")
    # globally_sharpened_map = prepare_sharpen_map(emmap_path,fsc_resolution=fsc_resolution,
    #                                        wilson_cutoff=wilson_cutoff, add_blur=add_blur,
    #                                        verbose=verbose,Cref=Cref)
    # 

    # Using the original emmap for refinement 
    
    #############################################################################
    # Stage 2b: Run servalcat to refine the reference model (either 
    # using the input PDB or the pseudo-atomic model)
    #############################################################################
    if verbose:
        print("."*80)
        print("Running model refinement\n")
    if skip_refine:
        if verbose: 
            tabbed_print.tprint("Skipping model refinements based on user input\n")
        refined_model_path = input_pdb_path
    else:
        if halfmap_paths is None:
            target_map = emmap_path
        else:
            target_map = halfmap_paths
        nyquist_resolution = 2*apix + 0.1
        refined_model_path = run_servalcat_iterative(model_path=input_pdb_path,  map_path=target_map,\
                    pseudomodel_refinement=pseudomodel_refinement, resolution=nyquist_resolution, num_iter=refmac_iter,\
                    refmac5_path=refmac5_path,verbose=verbose, hybrid_model_refinement=complete_model, \
                    final_chain_counts=final_chain_counts, cif_info=cif_info)
        
        if refined_model_path is None:
            tabbed_print.tprint("Problem running servalcat. Returning None")
            return None
        
    if os.path.exists(refined_model_path):
        bfactors = get_bfactors(refined_model_path)
            
        if verbose: 
            tabbed_print.tprint("ADP statistics for the refined model")
            tabbed_print.tprint("Mean B-factor: {}".format(np.mean(bfactors)))
            tabbed_print.tprint("Median B-factor: {}".format(np.median(bfactors)))
            tabbed_print.tprint("Max B-factor: {}".format(np.max(bfactors)))
            tabbed_print.tprint("Min B-factor: {}".format(np.min(bfactors)))
            ## If range of bfactors is too small then warn the user
            if max(bfactors)-min(bfactors) < 10:
                tabbed_print.tprint("Warning: The range of B-factors in the refined model is too small. Please check the model.")
                tabbed_print.tprint("Consider increasing the bfactor of the target map for refinement using the --add_blur option")
                tabbed_print.tprint("Current value used for add_blur = {}".format(add_blur))
        
        ## Now shift the refined bfactors to sharpen the emmap if required
        if not skip_refine:

            
            minimum_bfactor = 0    
            shifted_bfactors_structure, shift_value = shift_bfactors_by_probability(
                                        input_pdb=refined_model_path, probability_threshold=0.01, minimum_bfactor=minimum_bfactor)
            if verbose:
                tabbed_print.tprint("Shifting B-factor such that bfactor of p(<0.01) is {} (default)".format(minimum_bfactor))
                tabbed_print.tprint("Shifted B-factor by {}".format(shift_value))
            shifted_model_path = refined_model_path[:-4] + '_shifted_bfactors.pdb'
            shifted_bfactors_structure.write_pdb(shifted_model_path)

            if verbose:
                tabbed_print.tprint("Writing the shifted model to {}".format(shifted_model_path))
                # Print the statistics of the shifted model
                bfactors = get_bfactors(shifted_model_path)
                tabbed_print.tprint("ADP statistics for the shifted model")
                tabbed_print.tprint("Mean B-factor: {}".format(np.mean(bfactors)))
                tabbed_print.tprint("Median B-factor: {}".format(np.median(bfactors)))
                tabbed_print.tprint("Max B-factor: {}".format(np.max(bfactors)))
                tabbed_print.tprint("Min B-factor: {}".format(np.min(bfactors)))
        else:
            shifted_model_path = refined_model_path

    #############################################################################
    # Stage 3: Convert the refined model to a model-map using the 
    # run_refmap() function
    #############################################################################

    if verbose:
        print("."*80)
        print("Simulating model-map using refined structure factors\n")
    
    pseudomodel_modmap = run_refmap(model_path=shifted_model_path, emmap_path=emmap_path, mask_path=mask_path, verbose=verbose)
    
    #############################################################################
    # Stage 3a: If the user has specified symmetry, then apply the PG symmetry
    #############################################################################
    if pg_symmetry != "C1":
        if verbose:
            tabbed_print.tprint("Imposing a symmetry condition of {}".format(pg_symmetry))
        from locscale.include.symmetry_emda.symmetrize_map import symmetrize_map_emda
        from locscale.include.emmer.ndimage.map_utils import save_as_mrc
        sym = symmetrize_map_emda(emmap_path=pseudomodel_modmap,pg=pg_symmetry)
        symmetrised_modmap = pseudomodel_modmap[:-4]+"_{}_symmetry.mrc".format(pg_symmetry)
        save_as_mrc(map_data=sym, output_filename=symmetrised_modmap, apix=apix, origin=0, verbose=False)
        pseudomodel_modmap = symmetrised_modmap
    else:
        if verbose:
            tabbed_print.tprint("No symmetry condition imposed")
    
    #############################################################################
    # Stage 3b: If the user has specified a low pass filter cutoff then 
    # apply the low pass filter for model map
    #############################################################################
    if model_resolution is not None:
        if verbose:
            tabbed_print.tprint("Performing low pass filter on the Model Map with a cutoff: {} based on user input".format(model_resolution))
        from locscale.include.emmer.ndimage.filter import low_pass_filter
        from locscale.include.emmer.ndimage.map_utils import save_as_mrc
        
        pseudo_map_unfiltered_data = mrcfile.open(pseudomodel_modmap).data
        pseudo_map_filtered_data = low_pass_filter(im=pseudo_map_unfiltered_data, cutoff=model_resolution, apix=apix)
        
        filename = pseudomodel_modmap[:-4]+"_filtered.mrc"
        save_as_mrc(map_data=pseudo_map_filtered_data, output_filename=filename, apix=apix)
        
        pseudomodel_modmap = filename
    
    #############################################################################
    # Stage 4: Check and return the model-map
    #############################################################################

    # Collect pipeline intermediate files and output and dump them into a pickle file
    if verbose:
        print("."*80)
        print("Collecting intermediate files and dumping into a pickle file\n")
    preprocessing_pipeline_directory = os.path.dirname(emmap_path)
    if not complete_model:
        difference_mask_path = "not_used"
    intermediate_outputs = {
        "refined_model_path": refined_model_path,
        "shifted_model_path": shifted_model_path,
        "pseudomodel_modmap": pseudomodel_modmap,
        "mask_path": mask_path,
        "emmap_path": emmap_path,
        "input_pdb_path": input_pdb_path,
        "difference_mask_path": difference_mask_path,
        "preprocessing_pipeline_directory": preprocessing_pipeline_directory,
    }

    with open(os.path.join(preprocessing_pipeline_directory,"intermediate_outputs.pickle"), "wb") as f:
        pickle.dump(intermediate_outputs, f)
    
    if pseudomodel_modmap is None:
        tabbed_print.tprint("Problem simulating map from refined model. Returning None")
        return None
    else:
        tabbed_print.tprint("Successfully created model map")
        return pseudomodel_modmap
    


    
    
