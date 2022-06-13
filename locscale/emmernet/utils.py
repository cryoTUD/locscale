def download_emmernet_model_from_url(download_folder):
    import wget
   
    url_model_based_emmernet = "https://surfdrive.surf.nl/files/index.php/s/HxRLgoZFYQEbf8Z/download"
    wget.download(url_model_based_emmernet, download_folder)

def extract_tar_files_in_folder(tar_folder, use_same_folder=True):
    import tarfile
    import os
    if use_same_folder == 0:
        target_folder = tar_folder
    else:
        target_folder = os.path.dirname(tar_folder)

    for file in os.listdir(tar_folder):
        if file.endswith(".tar.gz"):
            print("Extracting: {}".format(file))
            tar = tarfile.open(os.path.join(tar_folder,file))
            tar.extractall(target_folder)
            tar.close()

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

def check_and_download_emmernet_model(verbose=False):
    ## Check if Emmernet model is downloaded
    import os
    import locscale

    emmernet_model_folder = os.path.join(os.path.dirname(locscale.__file__), "emmernet", "emmernet_models")
    path_exists = os.path.exists(emmernet_model_folder)
    MB_EMMERNET_MODEL_DOWNLOADED = os.path.exists(os.path.join(emmernet_model_folder, "EMmerNet_MBfa.hdf5"))
    MF_EMMERNET_MODEL_DOWNLOADED = os.path.exists(os.path.join(emmernet_model_folder, "EMmerNet_MFfa.hdf5"))
    ensemble_EMMERNET_MODEL_DOWNLOADED = os.path.exists(os.path.join(emmernet_model_folder, "EMmerNet_MBMF.hdf5"))

    emmernet_downloaded = path_exists and MB_EMMERNET_MODEL_DOWNLOADED and MF_EMMERNET_MODEL_DOWNLOADED and ensemble_EMMERNET_MODEL_DOWNLOADED

    if not emmernet_downloaded:
        if verbose:
            print("Emmernet model folder does not exist. Downloading model...")
        os.makedirs(emmernet_model_folder, exist_ok=True)
        download_emmernet_model_from_url(emmernet_model_folder)
        if verbose:
            print("Model downloaded")
        extract_tar_files_in_folder(emmernet_model_folder)
        if verbose:
            print("Model extracted")
    else:
        if verbose:
            print("Emmernet model folder exists: {}".format(emmernet_model_folder))
    
    return emmernet_model_folder

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
    output_emmap_folder = parsed_inputs["emmap_folder"]
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


    