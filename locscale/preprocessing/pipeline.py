
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
    from locscale.preprocessing.headers import run_pam, run_refmap, run_servalcat_iterative
    from locscale.preprocessing.prediction import predict_model_map_from_input_map
    from locscale.include.emmer.ndimage.map_utils import measure_mask_parameters, average_voxel_size
    from locscale.include.emmer.pdb.pdb_tools import find_wilson_cutoff, add_pseudoatoms_to_input_pdb
    from locscale.include.emmer.pdb.pdb_utils import get_bfactors, shift_bfactors_by_probability
    from locscale.utils.file_tools import RedirectStdoutToLogger, print_ADP_statistics, print_downward_arrow
    from locscale.utils.plot_tools import tab_print
    import mrcfile
    import pickle
    import numpy as np
    import os
    
    ###########################################################################
    # Extract the inputs from the dictionary
    ###########################################################################
    tabbed_print = tab_print(2)
    
    emmap_path = modmap_args['xyz_emmap_path']
    halfmap_paths = modmap_args['halfmap_paths']
    mask_path = modmap_args['mask_path_raw']
    pdb_path = modmap_args['model_coordinates']
    pseudomodel_method = modmap_args['pseudomodel_method']
    distance = modmap_args['distance']
    total_iterations = modmap_args['total_iterations']
    ref_resolution = modmap_args['ref_resolution']
    refmac_iterations = modmap_args['refmac_iterations']
    add_blur = modmap_args['add_blur']
    skip_refine = modmap_args['skip_refine']
    refmac5_path = modmap_args['refmac5_path']
    symmetry = modmap_args['symmetry']
    model_resolution = modmap_args['model_resolution']
    molecular_weight = modmap_args['molecular_weight']
    build_ca_only = modmap_args['build_ca_only']
    verbose = modmap_args['verbose']
    #Cref = modmap_args['Cref']
    complete_model = modmap_args['complete_model']
    averaging_window = modmap_args['averaging_window']
    mask_threshold = modmap_args['mask_threshold']
    cif_info = modmap_args['cif_info']
    modality = modmap_args['modality']

    if verbose:
        print("."*80)
        print("Running model-map generation pipeline \n")

    
    # if verbose:
        # print("Model map arguments: \n")
        # ## Print keys and values of dictionary in a nice format
        # for key, value in modmap_args.items():
        #     # if a value is numpy array then print its shape
        #     if isinstance(value, np.ndarray):
        #         value = value.shape
        #     if key == "Cref":
        #         # Print Cref shape
        #         if value is not None:
        #             print("{} : {}".format(key, value.shape))
        #         else:
        #             print("{} : {}".format(key, value))
        #     else:
        #         print("{:<20} : {}".format(key, value))

    #########################################################################
    # Open data files and collect required inputs
    # #######################################################################    
    emmap_mrc = mrcfile.open(emmap_path)
    apix = average_voxel_size(emmap_mrc.voxel_size)
    
    pam_bond_length = distance
    pam_method = pseudomodel_method
    total_iterations = total_iterations
    resolution = ref_resolution
    verbose = verbose
    ###########################################################################
    # Stage: Check the required number of atoms for the pseudomodel
    ###########################################################################
    if modality in ["pseudo_model_build_and_refine", "partial_model_input_build_and_refine"]:
        if molecular_weight is None:
            with RedirectStdoutToLogger(modmap_args['logger'], wait_message="Measuring input mask parameters"):
                num_atoms,mask_dims = measure_mask_parameters(mask_path,verbose=True, edge_threshold=mask_threshold)
            num_atoms = int(round(num_atoms / 1.55)) # adding correction factor based on comparison with atomic models within modelled regions
        else:
            avg_mass_per_atom = 13.14  #amu
            num_atoms = int(molecular_weight * 1000.0 / avg_mass_per_atom)
    ###########################################################################
    # Stage: Check if the user requires to build only Ca atoms
    ###########################################################################
    if build_ca_only:
        num_atoms = int(num_atoms/9)  ## Assuming 9 atoms per residue
        pam_bond_length = 3.8  ## Ca atom distances for secondary structures
        pam_method = 'gradient'  ## use this exclusively for Gradient
        if pam_method != 'gradient':
            print("Using gradient method for building pseudo-atomic model!\
                Not using user input:\t {}".format(pam_method))
    ###########################################################################
    # Stage : If user has not provided a PDB path then build a 
    # pseudomodel using the run_pam() routine else use the PDB path directly
    ###########################################################################
    if modality == "pseudo_model_build_and_refine":
        if verbose:
            print_statement = "a) Running pseudo-atomic model generator with {} atoms".format(num_atoms)
            print(print_statement)
            modmap_args['logger'].info(print_statement)
        input_pdb_path = run_pam(emmap_path=emmap_path, mask_path=mask_path, threshold=mask_threshold, num_atoms=num_atoms, 
                                method=pam_method, bl=pam_bond_length,total_iterations=total_iterations,verbose=True)
        pseudomodel_refinement = True
        modmap_args['logger'].info("Pseudo-atomic model generated at {}".format(input_pdb_path))
        if input_pdb_path is None:
            print("Problem running pseudo-atomic model generator. Returning None")
            modmap_args['logger'].error("Problem running pseudo-atomic model generator. Returning None")
            return None
        final_chain_counts = None
    elif modality == "partial_model_input_build_and_refine":
        pseudomodel_refinement = False
        if verbose:
            print_statement = "a) Running pseudo-atomic model generator to complete the user-provided PDB"
            print(print_statement)
            modmap_args['logger'].info(print_statement)
        integrated_structure, final_chain_counts, difference_mask_path = add_pseudoatoms_to_input_pdb(
            pdb_path=pdb_path, mask_path=mask_path, emmap_path=emmap_path,\
            averaging_window=averaging_window, pseudomodel_method=pam_method, pseudomodel_iteration=total_iterations, \
            mask_threshold=mask_threshold, fsc_resolution=ref_resolution, \
            return_chain_counts=True, return_difference_mask=True) 

        input_pdb_path = pdb_path[:-4] + '_integrated_pseudoatoms.cif'
        integrated_structure.make_mmcif_document().write_file(input_pdb_path)
        modmap_args['logger'].info("Integrated structure written to: {}".format(input_pdb_path))
        if not os.path.exists(input_pdb_path):
            print("Problem running pseudo-atomic model generator. Returning None")
            modmap_args['logger'].error("Problem running pseudo-atomic model generator. Returning None")
            return None
    elif modality == "full_model_input_refine_and_map" or modality == "full_model_input_no_refine":
        pseudomodel_refinement = False
        final_chain_counts = None
        if verbose:
            print_statement = "a) Using user-provided PDB for refinement"
            print(print_statement)
            modmap_args['logger'].info(print_statement)
        input_pdb_path = pdb_path
    elif modality == "predict_model_map":
        if verbose:
            print_statement = "a) Predicting model map from input map"
            print(print_statement)
            modmap_args['logger'].info(print_statement)
        pseudomodel_modmap = predict_model_map_from_input_map(modmap_args)
        return pseudomodel_modmap
    else:
        raise ValueError("Unknown modality: {}".format(modality))
        
        
    ###########################################################################
    # Stage: Refine the reference model usign servalcat
    ###########################################################################
            
    wilson_cutoff = find_wilson_cutoff(mask_path=mask_path, return_as_frequency=False, verbose=False)
    
    #############################################################################
    # Stage: Run servalcat to refine the reference model (either 
    # using the input PDB or the pseudo-atomic model)
    #############################################################################
    if verbose:
        print_downward_arrow(2)
        print_statement = "b) Running model refinement"
        print(print_statement)
        modmap_args['logger'].info(print_statement)
    
    
    if skip_refine:
        if verbose: 
            print_statement = "b) Skipping model refinements based on user input"
            print(print_statement)
            modmap_args['logger'].info(print_statement)
        refined_model_path = input_pdb_path
    else:
        if halfmap_paths is None:
            target_map = emmap_path
        else:
            target_map = halfmap_paths
        nyquist_resolution = 2*apix + 0.1
        refined_model_path = run_servalcat_iterative(model_path=input_pdb_path,  map_path=target_map,\
                    pseudomodel_refinement=pseudomodel_refinement, resolution=nyquist_resolution, num_iter=refmac_iterations,\
                    refmac5_path=refmac5_path,verbose=verbose, hybrid_model_refinement=complete_model, \
                    final_chain_counts=final_chain_counts, cif_info=cif_info)
    
        modmap_args['logger'].info("Refined model written to: {}".format(refined_model_path))
        if refined_model_path is None:
            modmap_args['logger'].error("Problem running servalcat. Returning None")
            print("Problem running servalcat. Returning None")
            return None
    
    if os.path.exists(refined_model_path):
        bfactors = get_bfactors(refined_model_path)
            
        with RedirectStdoutToLogger(modmap_args['logger'], wait_message="Calculating ADP statistics"):
            print_ADP_statistics(bfactors)
        ## If range of bfactors is too small then warn the user
        if max(bfactors)-min(bfactors) < 10:
            warn_message = "Warning: The range of B-factors in the refined model is too small. Please check the model."
            modmap_args['logger'].warning(warn_message)
            print(warn_message)
            
        
        ## Now shift the refined bfactors to sharpen the emmap if required
        if not skip_refine:
            minimum_bfactor = 0    
            with RedirectStdoutToLogger(modmap_args['logger'], wait_message="Shifting B-factors"):
                print("Shifting B-factor such that bfactor of p(<0.01) is {} (default)".format(minimum_bfactor))
                shifted_bfactors_structure, shift_value = shift_bfactors_by_probability(
                                            input_pdb=refined_model_path, probability_threshold=0.01, minimum_bfactor=minimum_bfactor)
                print("Shifted B-factor by {}".format(shift_value))
                
                extension = os.path.splitext(refined_model_path)[-1]
                shifted_model_path = refined_model_path.replace(extension, f"_shifted_bfactors{extension}")
                shifted_bfactors_structure.make_mmcif_document().write_file(shifted_model_path)
                print("Writing the shifted model to {}".format(shifted_model_path))
                # Print the statistics of the shifted model
                shifted_bfactors = get_bfactors(shifted_model_path)
                print_ADP_statistics(shifted_bfactors)
            if verbose:
                print_statement = "Shifted model written to: {}".format(shifted_model_path)
                print(print_statement)
                modmap_args['logger'].info(print_statement)
        else:
            shifted_model_path = refined_model_path

    #############################################################################
    # Stage: Convert the refined model to a model-map using run_refmap()
    #############################################################################

    if verbose:
        print_downward_arrow(2)
        print_statement = "c) Converting the refined model to a model-map"
        print(print_statement)
        modmap_args['logger'].info(print_statement)
    
    with RedirectStdoutToLogger(modmap_args['logger'], wait_message="Simulating model map"):
        pseudomodel_modmap = run_refmap(model_path=shifted_model_path, emmap_path=emmap_path, mask_path=mask_path, verbose=verbose)
    
    #############################################################################
    # Stage: If the user has specified symmetry, then apply the PG symmetry
    #############################################################################
    if symmetry != "C1":
        from locscale.include.symmetry_emda.symmetrize_map import symmetrize_map_emda
        from locscale.include.emmer.ndimage.map_utils import save_as_mrc
        
        if verbose:
            print_downward_arrow(2)
            print_statement = "d) Applying symmetry: {}".format(symmetry)
            print(print_statement)
            modmap_args['logger'].info(print_statement)
        
        with RedirectStdoutToLogger(modmap_args['logger'], wait_message="Applying symmetry"):
            sym = symmetrize_map_emda(emmap_path=pseudomodel_modmap,pg=symmetry)
            symmetrised_modmap = pseudomodel_modmap[:-4]+"_{}_symmetry.mrc".format(symmetry)
            save_as_mrc(map_data=sym, output_filename=symmetrised_modmap, apix=apix, origin=0, verbose=True)
            pseudomodel_modmap = symmetrised_modmap
    else:
        if verbose:
            print_statement = "d) No symmetry applied"
            print(print_statement)
            modmap_args['logger'].info(print_statement)
            
    
    #############################################################################
    # Stage: If the user has specified a low pass filter cutoff then 
    # apply the low pass filter for model map
    #############################################################################
    if model_resolution is not None:
        from locscale.include.emmer.ndimage.filter import low_pass_filter
        from locscale.include.emmer.ndimage.map_utils import save_as_mrc
        
        if verbose:
            print_downward_arrow(2)
            print_statement = "e) Applying low pass filter to the model map with a cutoff: {}".format(model_resolution)
            print(print_statement)
            modmap_args['logger'].info(print_statement)
        
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
        print_downward_arrow(2)
        print_statement = "Collecting pipeline intermediate files and output and dumping them into a pickle file"
        print(print_statement)
        modmap_args['logger'].info(print_statement)
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

    intermediate_outputs_pickle_path = os.path.join(preprocessing_pipeline_directory,"intermediate_outputs.pickle")
    with open(intermediate_outputs_pickle_path, "wb") as f:
        pickle.dump(intermediate_outputs, f)
    
    print_statement = f"-->{intermediate_outputs_pickle_path}<--"
    modmap_args['logger'].info(print_statement)
    print(print_statement)
    if pseudomodel_modmap is None:
        print("Problem simulating map from refined model. Returning None")
        modmap_args['logger'].error("Problem simulating map from refined model. Returning None")
        return None
    else:
        print("Successfully created model map")
        modmap_args['logger'].info("Successfully created model map")
        return pseudomodel_modmap
    


    
    
