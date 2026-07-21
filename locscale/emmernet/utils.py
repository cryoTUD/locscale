import numpy as np

def download_emmernet_model_from_url(download_folder):
    import os
    import urllib.request
    from tqdm import tqdm

    class DownloadProgressBar(tqdm):
        """Adapts tqdm to urlretrieve's reporthook(block_num, block_size, total_size)."""
        def update_to(self, block_num=1, block_size=1, total_size=None):
            if total_size is not None:
                self.total = total_size
            self.update(block_num * block_size - self.n)

    #url_emmernet_models = "https://zenodo.org/record/8211668/files/emmernet.tar.gz?download=1"
    url_emmernet_models = "https://surfdrive.surf.nl/files/index.php/s/txtZ7WiyRwixcrB/download"
    # urlretrieve needs a full path, not a directory: name the file explicitly rather
    # than letting it be inferred from a URL that ends in "?download=1".
    destination = os.path.join(download_folder, "emmernet.tar.gz")

    with DownloadProgressBar(unit="B", unit_scale=True, unit_divisor=1024,
                             miniters=1, desc="Downloading EMmerNet models") as progress:
        urllib.request.urlretrieve(url_emmernet_models, destination, reporthook=progress.update_to)

    return destination


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
    ax1.set_xlabel(r"Spatial frequency 1/$\AA$")
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
    
    # set an attribute to inform this is a feature_enhance run
    setattr(args, "run_type", "feature_enhance")


def check_emmernet_dependencies(verbose=False):
    try:
        import numpy as np
        import mrcfile
        import torch 
        import locscale
        
        if verbose:
            print("Emmernet dependencies are present")
    except ImportError: 
        raise 

def check_and_download_emmernet_model(verbose=False):
    ## Check if Emmernet model is downloaded
    import os
    import locscale

    emmernet_model_folder = os.path.join(os.path.dirname(locscale.__file__), "emmernet", "emmernet_models", "emmernet")
    path_exists = os.path.exists(emmernet_model_folder)
    EMMERNET_HIGH_CONTEXT_MODEL_DOWNLOADED = os.path.exists(os.path.join(emmernet_model_folder, "emmernet_highcontext.pt"))
    EMMERNET_LOW_CONTEXT_MODEL_DOWNLOADED = os.path.exists(os.path.join(emmernet_model_folder, "emmernet_lowcontext.pt"))
    

    emmernet_downloaded = path_exists and EMMERNET_HIGH_CONTEXT_MODEL_DOWNLOADED and EMMERNET_LOW_CONTEXT_MODEL_DOWNLOADED

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
    from locscale.emmernet.emmernet_functions import calibrate_variance
    
    input_emmap_path = parsed_inputs["emmap_path"]
    input_emmap_folder = os.path.dirname(os.path.dirname(input_emmap_path))
    processing_files_folder = emmernet_output["output_processing_files"]
    output_emmap_filename = parsed_inputs["outfile"]
    verbose = parsed_inputs["verbose"]
    monte_carlo = parsed_inputs["monte_carlo"]
    emmap, apix = load_map(input_emmap_path)
    output_emmap_folder = os.path.dirname(input_emmap_path)

    if monte_carlo:
        emmernet_output_mean = emmernet_output["output_predicted_map_mean"]
        emmernet_output_var = emmernet_output["output_predicted_map_var"]        
        #emmernet_output_var_calibrated = calibrate_variance(emmernet_output_var)
            
        assert emmap.shape == emmernet_output_mean.shape, "Emmernet output mean map shape does not match input map shape"
        assert emmap.shape == emmernet_output_var.shape, "Emmernet output var map shape does not match input map shape"
    else:
        emmernet_output_map = emmernet_output["output_predicted_map_mean"]
        assert emmap.shape == emmernet_output_map.shape, "Emmernet output map shape does not match input map shape"

    if verbose:
        print("."*80)
        print("Saving Emmernet output to {}".format(output_emmap_filename))
        
    if monte_carlo:
        # check if output filename has extension
        output_has_extension = len(os.path.splitext(output_emmap_filename)) > 1
        if not output_has_extension:
            output_emmap_filename = output_emmap_filename + ".mrc"
        extension_output_filename = os.path.splitext(output_emmap_filename)[1]
        output_filename_mean = os.path.join(input_emmap_folder, output_emmap_filename)
        output_filename_var = os.path.join(processing_files_folder, output_emmap_filename.replace(extension_output_filename, "_variance"+extension_output_filename))
        #output_filename_var_calibrated = output_emmap_filename.replace(extension_output_filename, "_var_calibrated"+extension_output_filename)
        output_filename_for_locscale = os.path.join(input_emmap_folder, output_emmap_filename.replace(extension_output_filename, "_baseline"+extension_output_filename))
        save_as_mrc(emmernet_output_mean, output_filename_mean, apix, verbose=verbose)
        save_as_mrc(emmernet_output_var, output_filename_var, apix, verbose=verbose)
        #save_as_mrc(emmernet_output_var_calibrated, output_filename_var_calibrated, apix, verbose=verbose)
        emmernet_output["output_filename_mean"] = output_filename_mean
        emmernet_output["output_filename_var"] = output_filename_var
        #emmernet_output["output_filename_var_calibrated"] = output_filename_var_calibrated
        emmernet_output["output_filename_for_locscale"] = os.path.join(output_emmap_folder, output_filename_for_locscale)
        emmernet_output["reference_map_for_locscale"] = output_filename_mean
    else:
        save_as_mrc(emmernet_output_map, output_emmap_filename, apix, verbose=verbose)
        emmernet_output["output_filename"] = output_emmap_filename

        output_has_extension = len(os.path.splitext(output_emmap_filename)) > 1
        if not output_has_extension:
            output_emmap_filename = output_emmap_filename + ".mrc"
        extension_output_filename = os.path.splitext(output_emmap_filename)[1]
        # No Monte Carlo variance without dropout sampling: the single
        # prediction doubles as the reference map LocScale scales against.
        emmernet_output["output_filename_for_locscale"] = output_emmap_filename.replace(
            extension_output_filename, "_baseline" + extension_output_filename)
        emmernet_output["reference_map_for_locscale"] = output_emmap_filename

    return emmernet_output

def load_calibrator():
    from locscale.utils.file_tools import get_locscale_path
    import pickle 
    import os 
    
    locscale_path = get_locscale_path()
    regressor_path = os.path.join(locscale_path, "locscale", "utils", "calibrator_locscale_target_seed_42.pickle")
    calibrator = pickle.load(open(regressor_path, "rb"))
    
    return calibrator
    
def symmetrise_if_needed(input_dictionary, output_dictionary,):
    import os
    symmetry = input_dictionary["symmetry"]
    twist = input_dictionary.get("twist")
    rise = input_dictionary.get("rise")
    n_steps = input_dictionary.get("n_steps")

    if symmetry != "C1" or twist is not None:
        import torch 
        from locscale.include.symmetry_emda.symmetrize_map import symmetrize_map_emda
        from locscale.include.emmer.ndimage.map_utils import save_as_mrc, load_map
        from locscale.utils.file_tools import RedirectStdoutToLogger

        input_map_shape = input_dictionary["input_map_shape"][0]

        if input_map_shape < 350 and torch.cuda.is_available():
            device = torch.device("cuda")
        else:
            device = torch.device("cpu")


        verbose = input_dictionary["verbose"]
        map_to_symmetrise = output_dictionary["output_predicted_map_mean"]
        processing_files_folder = input_dictionary["output_processing_files"]
        # save the non-symmetrised map
        apix = input_dictionary["apix"]
        unsymmetrised_map_path = os.path.join(processing_files_folder, "unsymmetrised_mean_map.mrc")
        save_as_mrc(map_data=map_to_symmetrise, output_filename=unsymmetrised_map_path, apix=apix, origin=0, verbose=True)

        _, apix = load_map(unsymmetrised_map_path)
        if verbose:
            if twist is not None:
                print_statement = "Applying helical symmetry: twist={}, rise={}, point group={}".format(twist, rise, symmetry)
            else:
                print_statement = "Applying symmetry: {}".format(symmetry)
            print(print_statement)
            input_dictionary['logger'].info(print_statement)

        with RedirectStdoutToLogger(input_dictionary['logger'], wait_message="Applying symmetry"):
            sym = symmetrize_map_emda(emmap_path=unsymmetrised_map_path, pg=symmetry, twist=twist, rise=rise, n_steps=n_steps, device=device)
            symmetrised_map = unsymmetrised_map_path[:-4]+"_{}_symmetry.mrc".format(symmetry)
            save_as_mrc(map_data=sym, output_filename=symmetrised_map, apix=apix, origin=0, verbose=True)

        output_dictionary["output_predicted_map_mean_non_symmetrised"] = map_to_symmetrise
        output_dictionary["output_predicted_map_mean"] = sym

        return output_dictionary
    else:
        return output_dictionary
            
def compute_calibrated_probabilities(locscale_path, mean_prediction_path, variance_prediction_path, mask_path, n_samples=15):
    from locscale.emmernet.emmernet_functions import load_smoothened_mask
    from locscale.emmernet.utils import load_calibrator
    from locscale.include.emmer.ndimage.map_utils import load_map
    import numpy as np
    import os 
    
    #http://www.ltcconline.net/greenl/courses/201/estimation/smallConfLevelTable.htm
    z_target_for_each_CI = {70 : 1.04, 80: 1.28, 90: 1.645, 95: 1.96, 99: 2.58}
    
    
    locscale_map, apix = load_map(locscale_path)
    mean_prediction, _ = load_map(mean_prediction_path)
    variance_prediction, _ = load_map(variance_prediction_path)
    variance_mask, _ = load_map(mask_path)
    variance_mask = variance_mask > 0.5
        
    locscale_masked = locscale_map[variance_mask]
    mean_masked = mean_prediction[variance_mask]
    variance_masked = variance_prediction[variance_mask]
    
    standard_deviation_masked = np.sqrt(variance_masked)
    standard_error_masked = standard_deviation_masked / np.sqrt(n_samples)
    
    calibrator = load_calibrator()
    calibrated_standard_error = calibrator.predict(standard_error_masked)
    
    # compute the z-scores
    z_scores = (locscale_masked - mean_masked) / calibrated_standard_error
    
    # compute the probabilities for different confidence intervals
    observed_probabilities = {}
    for ci in z_target_for_each_CI:
        z_target = z_target_for_each_CI[ci]
        observed_probability = np.sum(np.abs(z_scores) < z_target)
        observed_probabilities[ci] = observed_probability / len(z_scores)
        
    
    return observed_probabilities


def plot_binned_correlation(xarray, yarray, num_bins=50, ci = 0.95, figsize_cm=(8, 8), plot_diagonal=True):
    import matplotlib.pyplot as plt
    import scipy.stats as st
    import numpy as np
                
    figsize = (figsize_cm[0] / 2.54, figsize_cm[1] / 2.54)
    fig, ax = plt.subplots(figsize=figsize)

    import warnings
    warnings.filterwarnings('ignore')
    xarray = np.array(xarray)
    yarray = np.array(yarray)
    
    # Binning by the x axis data
    bins = np.linspace(xarray.min(), xarray.max(), num_bins)
    bin_indices = np.digitize(xarray, bins)

    # Compute the statistics for each bin
    bin_not_empty = lambda i: len(xarray[bin_indices == i]) > 0
    
    xarray_bin_means = [xarray[bin_indices == i].mean() for i in range(len(bins)) if bin_not_empty(i)]
    yarray_bin_means = [yarray[bin_indices == i].mean() for i in range(len(bins)) if bin_not_empty(i)]
    yarray_bin_stds = [yarray[bin_indices == i].std() for i in range(len(bins)) if bin_not_empty(i)]
    yarray_bin_nums = [len(yarray[bin_indices == i]) for i in range(len(bins)) if bin_not_empty(i)]
    z_score = st.norm.ppf(ci)
    yarray_standard_errors = [z_score * yarray_bin_stds[i] / np.sqrt(yarray_bin_nums[i]) for i in range(len(xarray_bin_means))]
    
    # convert to numpy arrays
    xarray_bin_means = np.array(xarray_bin_means)
    yarray_bin_means = np.array(yarray_bin_means)
    yarray_standard_errors = np.array(yarray_standard_errors)
    
    # Find max yarray and min yarray
    yarray_top = yarray_bin_means + yarray_standard_errors
    yarray_bottom = yarray_bin_means - yarray_standard_errors
        
    ax.plot(xarray_bin_means, yarray_bin_means, color='blue', marker='o')

    max_x, max_y = np.max(xarray_bin_means), np.max(yarray_bin_means)
    min_x, min_y = np.min(xarray_bin_means), np.min(yarray_bin_means)
    min_val, max_val = np.min([min_x, min_y]), np.max([max_x, max_y])
    if plot_diagonal:
        # plot diagonal line for reference
        ax.plot([min_val, max_val], [min_val, max_val], color='red', linestyle='--')
    # shade the area between the standard errors
    ax.fill_between(xarray_bin_means, yarray_bottom, yarray_top, color='skyblue', alpha=0.4, label=f'ci:{ci}')
    
    plt.tight_layout()
    
    return fig, ax

def compute_reliability_curve(locscale_path, mean_prediction_path, variance_prediction_path, mask_path, n_samples=15):
    from locscale.emmernet.utils import load_calibrator
    from locscale.include.emmer.ndimage.map_utils import load_map
    import numpy as np
    import os     
    
    locscale_map, apix = load_map(locscale_path)
    mean_prediction, _ = load_map(mean_prediction_path)
    variance_prediction, _ = load_map(variance_prediction_path)
    variance_mask, _ = load_map(mask_path)
    variance_mask = variance_mask > 0.5
    
    locscale_masked = locscale_map[variance_mask]
    mean_masked = mean_prediction[variance_mask]
    variance_masked = variance_prediction[variance_mask]
    
    standard_deviation_masked = np.sqrt(variance_masked)
    standard_error_masked = standard_deviation_masked / np.sqrt(n_samples)
    
    calibrator = load_calibrator()
    calibrated_standard_error = calibrator.predict(standard_error_masked)
    
    absolute_residual = np.abs(locscale_masked - mean_masked)
    
    fig, ax = plot_binned_correlation(calibrated_standard_error, absolute_residual, num_bins=128, ci=0.95, figsize_cm=(8, 8))
    # modify plot 
    ax.set_xlabel("Calibrated Standard Error")
    ax.set_ylabel("Absolute Residual")
    ax.set_title("Reliability Curve")
    
    xmin, xmax = ax.get_xlim()
    ymin, ymax = ax.get_ylim()
    maxval, minval = np.max([xmax, ymax]), np.min([xmin, ymin])
    
    ax.set_xlim(minval, maxval)
    ax.set_ylim(minval, maxval)
    
    return fig, ax 