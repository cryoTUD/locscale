import numpy as np

def download_emmernet_model_from_url(download_folder):
    import wget
   
    #url_model_based_emmernet = "https://surfdrive.surf.nl/files/index.php/s/HxRLgoZFYQEbf8Z/download"    # OLD SURFdrive link
    url_model_based_emmernet = "https://zenodo.org/record/6651995/files/emmernet.tar.gz?download=1"    # NEW Zenodo link
    wget.download(url_model_based_emmernet, download_folder)

def extract_tar_files_in_folder(tar_folder, use_same_folder=True):
    import tarfile
    import os
    if use_same_folder:
        target_folder = tar_folder
    else:
        target_folder = os.path.dirname(tar_folder)

    for file in os.listdir(tar_folder):
        if file.endswith(".tar.gz"):
            print("\nExtracting: {}".format(file))
            tar = tarfile.open(os.path.join(tar_folder,file))
            tar.extractall(target_folder)
            tar.close()

def compute_local_phase_correlations(target_cubes, predicted_cubes, apix, temp_folder=None):
    import os
    from tqdm import tqdm
    from locscale.include.emmer.ndimage.fsc_util import calculate_phase_correlation_maps
    from locscale.include.emmer.ndimage.profile_tools import frequency_array
    phase_correlations_all = []
    for i in tqdm(range(len(predicted_cubes)), desc="Calculating phase correlations"):
        predicted_cube = predicted_cubes[i]
        target_cube = target_cubes[i]
        if target_cube.sum() > 6000:
            phase_correlation_cube = calculate_phase_correlation_maps(predicted_cube, target_cube)
            phase_correlations_all.append(phase_correlation_cube[1:])
    phase_correlations_all = np.array(phase_correlations_all)
    freq = frequency_array(phase_correlation_cube, apix)
    # save the numpy arrays 
    if temp_folder is not None:
        np.save(os.path.join(temp_folder, "target_cubes.npy"), target_cubes)
        np.save(os.path.join(temp_folder, "predicted_cubes.npy"), predicted_cubes)
        np.save(os.path.join(temp_folder, "phase_correlations_all.npy"), phase_correlations_all)
        np.save(os.path.join(temp_folder, "freq.npy"), freq)
    
    return phase_correlations_all, freq

def plot_phase_correlations(phase_correlations_all, freq):
    import matplotlib.pyplot as plt
    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    ax1.grid(False)
    ax2 = ax1.twiny()
    for phase_correlation in phase_correlations_all:
        ax1.plot(freq[1:], phase_correlation, color="black", alpha=0.1)
    # Plot the mean
    ax1.plot(freq[1:], phase_correlations_all.mean(axis=0), color="red", linewidth=2)
    ax1.set_xlabel("Spatial frequency 1/$\AA$")
    ax1.set_ylabel("Phase correlation")
    ax2.set_xticks(ax1.get_xticks())
    ax2.set_xbound(ax1.get_xbound())
    ax2.set_xticklabels([round(1/x,1) for x in ax1.get_xticks()])
    ax2.set_xlabel(r'Resolution $(\AA)$')
    plt.ylim(-0.5,1.2)
    # add Y tick labels as [0, 0.5, 1]
    ax1.set_yticks([0, 0.5, 1])
    plt.tight_layout()
    return fig


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
            print("\nEmmernet model folder does not exist. Downloading model...\n")
        os.makedirs(emmernet_model_folder, exist_ok=True)
        download_emmernet_model_from_url(emmernet_model_folder)
        if verbose:
            print("Model downloaded\n")
        extract_tar_files_in_folder(emmernet_model_folder, use_same_folder=True)
        if verbose:
            print("Model extracted\n")
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
    verbose = parsed_inputs["verbose"]
    monte_carlo = parsed_inputs["monte_carlo"]
    
    emmap, apix = load_map(input_emmap_path)

    if monte_carlo:
        emmernet_output_mean = emmernet_output["output_mean"]
        emmernet_output_var = emmernet_output["output_var"]
        
        assert emmap.shape == emmernet_output_mean.shape, "Emmernet output mean map shape does not match input map shape"
        assert emmap.shape == emmernet_output_var.shape, "Emmernet output var map shape does not match input map shape"
    else:
        emmernet_output_map = emmernet_output["output"]
        #emmernet_output_charge_density = emmernet_output["output_charge_density"]
        assert emmap.shape == emmernet_output_map.shape, "Emmernet output map shape does not match input map shape"

    if verbose:
        print("."*80)
        print("Saving Emmernet output to {}".format(output_emmap_filename))
        
    if monte_carlo:
        extension_output_filename = os.path.splitext(output_emmap_filename)[1]
        output_filename_mean = output_emmap_filename.replace(extension_output_filename, "_mean"+extension_output_filename)
        output_filename_var = output_emmap_filename.replace(extension_output_filename, "_var"+extension_output_filename)
        save_as_mrc(emmernet_output_mean, output_filename_mean, apix, verbose=verbose)
        save_as_mrc(emmernet_output_var, output_filename_var, apix, verbose=verbose)
    else:
        save_as_mrc(emmernet_output_map, output_emmap_filename, apix, verbose=verbose)
        #extension_output_filename = os.path.splitext(output_emmap_filename)[1]
        #charge_density_filename = output_emmap_filename.replace(extension_output_filename, "_charge_density"+extension_output_filename)
        #save_as_mrc(emmernet_output_charge_density, charge_density_filename, apix, verbose=verbose)


    