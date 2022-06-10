
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
    from locscale.preprocessing.headers import run_FDR, run_pam, run_refmac_servalcat, run_refmap, prepare_sharpen_map, is_pseudomodel
    from locscale.include.emmer.ndimage.map_utils import measure_mask_parameters, average_voxel_size
    from locscale.include.emmer.pdb.pdb_tools import find_wilson_cutoff
    from locscale.utils.plot_tools import tab_print
    import mrcfile
    
    tabbed_print = tab_print(2)
    
    emmap_path = modmap_args['emmap_path']
    mask_path = modmap_args['mask_path']
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

    if verbose:
        print("."*80)
        print("Running model-map generation pipeline \n")


    if verbose:
        tabbed_print.tprint("Model map arguments: \n")
        ## Print keys and values of dictionary in a nice format
        for key, value in modmap_args.items():
            tabbed_print.tprint("{:<20} : {}".format(key, value))

    ################### START PIPELINE ###################
    
    emmap_mrc = mrcfile.open(emmap_path)
    apix = average_voxel_size(emmap_mrc.voxel_size)
    pam_bond_length = pam_distance
    pam_method = pseudomodel_method
    pam_iteration = pam_iteration
    
    resolution = fsc_resolution
    verbose = verbose
    
    if molecular_weight is None:
        num_atoms,mask_dims = measure_mask_parameters(mask_path,verbose=False)
    else:
        avg_mass_per_atom = 13.14  #amu
        num_atoms = int(molecular_weight * 1000.0 / avg_mass_per_atom)

    if build_ca_only:
        num_atoms = int(num_atoms/9)  ## Assuming 9 atoms per residue
        pam_bond_length = 3.8  ## Ca atom distances for secondary structures
        pam_method = 'gradient'  ## use this exclusively for Gradient
        if pam_method != 'gradient':
            tabbed_print.tprint("Using gradient method for building pseudo-atomic model!\
                Not using user input:\t {}".format(pam_method))
    
    if pdb_path is None:
        if verbose:
            print("."*80)
            print("You have not entered a PDB path, running pseudo-atomic model generator!")
        input_pdb_path = run_pam(emmap_path=emmap_path, mask_path=mask_path, threshold=1, num_atoms=num_atoms, 
                                   method=pam_method, bl=pam_bond_length,total_iterations=pam_iteration,verbose=verbose)
        if input_pdb_path is None:
            print("Problem running pseudo-atomic model generator. Returning None")
            return None
    else:
        if verbose:
            print("."*80)
            print("Using user-provided PDB path: {}".format(pdb_path))
        input_pdb_path = pdb_path
    
    if is_pseudomodel(input_pdb_path):
        only_bfactor_refinement = True
    else:
        only_bfactor_refinement = False
            
    wilson_cutoff = find_wilson_cutoff(mask_path=mask_path, return_as_frequency=False, verbose=False)
    
    if verbose:
        print("."*80)
        print("Preparing target map for refinement\n")
    globally_sharpened_map = prepare_sharpen_map(emmap_path,fsc_resolution=fsc_resolution,
                                            wilson_cutoff=wilson_cutoff, add_blur=add_blur, verbose=verbose)
    
    if verbose:
        print("."*80)
        print("Running Refmac refinement\n")
    if skip_refine:
        if verbose: 
            tabbed_print.tprint("Skipping REFMAC refinements based on user input\n")
        refined_model_path = input_pdb_path
    else:
        refined_model_path = run_refmac_servalcat(model_path=input_pdb_path,  map_path=globally_sharpened_map,\
                    only_bfactor_refinement=only_bfactor_refinement, resolution=resolution, num_iter=refmac_iter,
                    refmac5_path=refmac5_path,verbose=verbose)
        if refined_model_path is None:
            tabbed_print.tprint("Problem running REFMAC. Returning None")
            return None
        
    if verbose:
        print("."*80)
        print("Simulating model-map using refined structure factors\n")
    pseudomodel_modmap = run_refmap(model_path=refined_model_path, emmap_path=emmap_path, mask_path=mask_path, verbose=verbose)
    
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
    
    
    
    if pseudomodel_modmap is None:
        tabbed_print.tprint("Problem simulating map from refined model. Returning None")
        return None
    else:
        tabbed_print.tprint("Successfully created model map")
        return pseudomodel_modmap
    


    
    
