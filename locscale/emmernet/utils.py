def check_emmernet_inputs(args):
    '''
    Check user inputs for errors and conflicts

    Parameters
    ----------
    args : TYPE
        DESCRIPTION.

    Returns
    -------
    None.

    '''
    
    import mrcfile
    import os
    from textwrap import fill

    ## Check input files
    emmap_absent = True
    if args.emmap_path is not None:
        if os.path.exists(args.emmap_path):
            emmap_absent = False
    
    half_maps_absent = True
    if args.halfmap_paths is not None:
        halfmap1_path = args.halfmap_paths[0]
        halfmap2_path = args.halfmap_paths[1]
        if os.path.exists(halfmap1_path) and os.path.exists(halfmap2_path):
            half_maps_absent = False
    
    
    if args.outfile is None:
        print(fill("You have not entered a filename for EMmerNet output. Using a standard output file name: emmernet_prediction.mrc. \
            Any file with the same name in the current directory will be overwritten", 80))
        print("\n")

        outfile = [x for x in vars(args) if x=="outfile"]
        
        setattr(args, outfile[0], "emmernet_prediction.mrc")


def check_emmernet_dependencies(verbose=False):
    try:
        import numpy as np
        import mrcfile
        import tensorflow as tf
        import keras
        import locscale
        
        if verbose:
            print("Emmernet dependencies are present")
    except ImportError: 
        raise 

def check_and_save_output(parsed_inputs, emmernet_output):
    '''
    Check if the output file is present and save the output if it is not.

    Parameters
    ----------
    parsed_inputs : dictionary
        .
    emmernet_output : dictionary
        .

    Returns
    -------
    None.

    '''
    import os
    from locscale.include.emmer.ndimage.map_utils import save_as_mrc, load_map

    input_emmap_path = parsed_inputs["emmap_path"]
    input_emmap_folder = os.path.dirname(input_emmap_path)
    output_emmap_filename = parsed_inputs["outfile"]
    output_emmap_folder = os.path.dirname(output_emmap_filename)
    verbose = parsed_inputs["verbose"]

    if output_emmap_folder is not None:
        if not os.path.exists(output_emmap_folder):
            output_emmap_folder = input_emmap_folder
    
    output_emmap_path = os.path.join(output_emmap_folder, output_emmap_filename)
    
    emmap, apix = load_map(input_emmap_path)

    emmernet_output_map = emmernet_output["output"]

    assert emmap.shape == emmernet_output_map.shape, "Emmernet output map shape does not match input map shape"

    if verbose:
        print("."*80)
        print("Saving Emmernet output to {}".format(output_emmap_path))
        

    save_as_mrc(emmernet_output_map, output_emmap_path, apix, verbose=verbose)


    