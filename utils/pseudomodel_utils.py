
def get_modmap_from_pseudomodel(emmap_path, mask_path, pseudomodel_method, pam_distance, pam_iteration, fsc_resolution, verbose):
    '''
    Function to generate a model map using pseudo-atomic model

    Parameters
    ----------
    emmap_path : TYPE
        DESCRIPTION.
    mask_path : TYPE
        DESCRIPTION.
    pseudomodel_method : TYPE
        DESCRIPTION.
    pam_distance : TYPE
        DESCRIPTION.
    pam_iteration : TYPE
        DESCRIPTION.
    fsc_resolution : TYPE
        DESCRIPTION.
    verbose : TYPE
        DESCRIPTION.

    Returns
    -------
    pseudomodel_modmap : TYPE
        DESCRIPTION.

    '''
    from locscale_headers import run_FDR, measure_mask_parameters, run_pam, run_refmac, run_refmap, prepare_sharpen_map
    from emmer.pdb.pdb_tools import find_wilson_cutoff
    import mrcfile
    
    emmap_mrc = mrcfile.open(emmap_path)
       
    pam_bond_length = pam_distance
    pam_method = pseudomodel_method
    pam_iteration = pam_iteration
    
    resolution = fsc_resolution
    verbose = verbose
    
    
    num_atoms,mask_dims = measure_mask_parameters(mask_path,verbose=verbose)
    
    pseudomodel_path = run_pam(emmap_path=emmap_path, mask_path=mask_path, threshold=1, num_atoms=num_atoms, 
                               method=pam_method, bl=pam_bond_length,total_iterations=pam_iteration,verbose=verbose)
    if pseudomodel_path is None:
        print("Problem running pseudo-atomic model generator. Returning None")
        return None
    
    wilson_cutoff = find_wilson_cutoff(mask_path=mask_path, return_as_frequency=False)
    
    globally_sharpened_map = prepare_sharpen_map(emmap_path, wilson_cutoff=wilson_cutoff)
    
    refined_model_path = run_refmac(model_path=pseudomodel_path, model_name=pseudomodel_path[:-4], 
                                    map_path=globally_sharpened_map, resolution=resolution, maskdims=mask_dims,verbose=verbose)
    if refined_model_path is None:
        print("Problem running REFMAC. Returning None")
        return None
    
    #emmap_path, mask_path = run_mapmask(args.em_map), run_mapmask(mask_path)
    #pseudomodel_modmap,new_emmap_path,new_mask_path = run_refmap2(model_path=refined_model_path, 
                                                                 #emmap_path=args.em_map, 
                                                                 #mask_path=mask_path, verbose=verbose)
    
    pseudomodel_modmap = run_refmap(model_path=refined_model_path, emmap_path=emmap_path, mask_path=mask_path, verbose=verbose)
    
    
    
    if pseudomodel_modmap is None:
        print("Problem simulating map from refined model. Returning None")
        return None
    else:
        print("Successfully created model map")
        return pseudomodel_modmap
    


    
    