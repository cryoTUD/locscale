import numpy as np
import math

## MAP PROCESSING FUNCTIONS REQUIRED FOR LOCSCALE 

def shift_map_to_zero_origin(emmap_path):
    '''
    Determines the map origin from header file and changes it to zero

    Parameters
    ----------
    emmap_path : str
        DESCRIPTION.

    Returns
    -------
    shift_vector : numpy.ndarray (len=3)

    '''    
    import mrcfile
    from locscale.include.emmer.ndimage.map_utils import save_as_mrc
    
    target_origin = np.array([0,0,0])
    voxel_size = np.array(mrcfile.open(emmap_path).voxel_size.tolist())
    current_origin = np.array(mrcfile.open(emmap_path).header.origin.tolist()) 
    
    emmap_data = mrcfile.open(emmap_path).data
    
    output_file = emmap_path
    save_as_mrc(map_data=emmap_data, output_filename=emmap_path, apix=voxel_size, origin=0)
    
    shift_vector = target_origin - current_origin
    return shift_vector

def get_spherical_mask(emmap):
    from locscale.utils.general import pad_or_crop_volume
    
    mask = np.zeros(emmap.shape)

    if mask.shape[0] == mask.shape[1] and mask.shape[0] == mask.shape[2] and mask.shape[1] == mask.shape[2]:
        rad = mask.shape[0] // 2
        z,y,x = np.ogrid[-rad: rad+1, -rad: rad+1, -rad: rad+1]
        mask = (x**2+y**2+z**2 <= rad**2).astype(np.int_).astype(np.int8)
        mask = pad_or_crop_volume(mask,emmap.shape)
        mask = (mask==1).astype(np.int8)
    else:
        mask += 1
        mask = mask[0:mask.shape[0]-1, 0:mask.shape[1]-1, 0:mask.shape[2]-1]
        mask = pad_or_crop_volume(emmap, (emmap.shape), pad_value=0)
    
    return mask

def put_scaled_voxels_back_in_original_volume_including_padding(sharpened_vals, masked_indices, map_shape):
    map_scaled = np.zeros(np.prod(map_shape))
    map_scaled[masked_indices] = sharpened_vals
    map_scaled = map_scaled.reshape(map_shape)

    return map_scaled

def pad_or_crop_volume(vol, dim_pad=None, pad_value = None, crop_volume=False):
    from locscale.utils.math_tools import round_up_proper
    if (dim_pad == None):
        return vol
    else:
        dim_pad = np.round(np.array(dim_pad)).astype('int')
        #print(dim_pad)

        if pad_value == None:
            pad_value = 0

        if (dim_pad[0] <= vol.shape[0] or dim_pad[1] <= vol.shape[1] or dim_pad[2] <= vol.shape[2]):
            crop_volume = True

        if crop_volume:
            k_start = round_up_proper(vol.shape[0]/2-dim_pad[0]/2)
            k_end = round_up_proper(vol.shape[0]/2+dim_pad[0]/2)
            j_start = round_up_proper(vol.shape[1]/2-dim_pad[1]/2)
            j_end = round_up_proper(vol.shape[1]/2+dim_pad[1]/2)
            i_start = round_up_proper(vol.shape[2]/2-dim_pad[2]/2)
            i_end = round_up_proper(vol.shape[2]/2+dim_pad[2]/2)
            crop_vol = vol[k_start:k_end, :, :]
            crop_vol = crop_vol[:, j_start:j_end, :]
            crop_vol = crop_vol[:, :, i_start:i_end]

            return crop_vol

        else:
            k_start = round_up_proper(dim_pad[0]/2-vol.shape[0]/2)
            k_end = round_up_proper(dim_pad[0]/2-vol.shape[0]/2)
            j_start = round_up_proper(dim_pad[1]/2-vol.shape[1]/2)
            j_end = round_up_proper(dim_pad[1]/2-vol.shape[1]/2)
            i_start = round_up_proper(dim_pad[2]/2-vol.shape[2]/2)
            i_end = round_up_proper(dim_pad[2]/2-vol.shape[2]/2)
            
            pad_vol = np.pad(vol, ((k_start, k_end ), (0,0), (0,0) ), 'constant', constant_values=(pad_value,))
            pad_vol = np.pad(pad_vol, ((0,0), (j_start, j_end ), (0,0)), 'constant', constant_values=(pad_value,))
            pad_vol = np.pad(pad_vol, ((0,0), (0,0), (i_start, i_end )), 'constant', constant_values=(pad_value,))
            
            if pad_vol.shape[0] != dim_pad[0] or pad_vol.shape[1] != dim_pad[1] or pad_vol.shape[2] != dim_pad[2]:
                print("Requested pad volume shape {} not equal to the shape of the padded volume returned{}. Input map shape might be an odd sized map.".format(dim_pad, pad_vol.shape))
                

            return pad_vol

def compute_padding_average(vol, mask):
    mask = (mask == 1).astype(np.int8)
    #inverted_mask = np.logical_not(mask)
    average_padding_intensity = np.mean(np.ma.masked_array(vol, mask))
    return average_padding_intensity

def get_xyz_locs_and_indices_after_edge_cropping_and_masking(mask, wn):
    mask = np.copy(mask)
    nk, nj, ni = mask.shape

    kk, jj, ii = np.indices((mask.shape))
    kk_flat = kk.ravel()
    jj_flat = jj.ravel()
    ii_flat = ii.ravel()

    mask_bin = np.array(mask.ravel(), dtype=np.bool)
    indices = np.arange(mask.size)
    masked_indices = indices[mask_bin]
    cropped_indices = indices[(wn / 2 <= kk_flat) & (kk_flat < (nk - wn / 2)) &
                              (wn / 2 <= jj_flat) & (jj_flat < (nj - wn / 2)) &
                              (wn / 2 <= ii_flat) & (ii_flat < (ni - wn / 2))]

    cropp_n_mask_ind = np.intersect1d(masked_indices, cropped_indices)

    xyz_locs = np.column_stack((kk_flat[cropp_n_mask_ind], jj_flat[cropp_n_mask_ind], ii_flat[cropp_n_mask_ind]))

    return xyz_locs, cropp_n_mask_ind, mask.shape

def check_for_window_bleeding(mask,wn):
    from locscale.utils.general import get_xyz_locs_and_indices_after_edge_cropping_and_masking
    masked_xyz_locs, masked_indices, mask_shape = get_xyz_locs_and_indices_after_edge_cropping_and_masking(mask, 0)

    zs, ys, xs = masked_xyz_locs.T
    nk, nj, ni = mask_shape
    #print(xs.shape, ys.shape, zs.shape)
    #print(nk,nj,ni)
    #print(wn)

    if xs.min() < wn / 2 or xs.max() > (ni - wn / 2) or \
    ys.min() < wn / 2 or ys.max() > (nj - wn / 2) or \
    zs.min() < wn / 2 or zs.max() > (nk - wn / 2):
        window_bleed = True
    else:
        window_bleed = False

    return window_bleed

def normalise_intensity_levels(from_emmap, to_levels=[0,1]):
    normalise_between_zero_one = (from_emmap - from_emmap.min()) / (from_emmap.max() - from_emmap.min())
    to_levels = np.array(to_levels)
    
    min_value = to_levels.min()
    max_value = to_levels.max()
    scale_factor = max_value-min_value
    
    normalised = min_value + normalise_between_zero_one * scale_factor
    
    return normalised    

def shift_radial_profile(from_emmap, to_emmap):
    '''
    To shift the radial profile of one emmap so that DC power matches another emmap

    Parameters
    ----------
    from_emmap : numpy.ndarray
        DESCRIPTION.
    to_emmap : numpy.ndarray
        DESCRIPTION.

    Returns
    -------
    emmap_shifted_rp : numpy.ndarray

    '''
    from locscale.include.emmer.ndimage.profile_tools import compute_radial_profile, plot_radial_profile, frequency_array
    from locscale.include.emmer.ndimage.map_tools import set_radial_profile, compute_scale_factors
    
    rp_from_emmap = compute_radial_profile(from_emmap)
    rp_to_emmap, radii = compute_radial_profile(to_emmap, return_indices=True)
    
    DC_power_diff = rp_to_emmap.max() - rp_from_emmap.max()
    print((DC_power_diff))
    new_rp_from_emmap = rp_from_emmap * 20
    scale_factors = compute_scale_factors(rp_from_emmap, new_rp_from_emmap)
    emmap_shifted_rp = set_radial_profile(from_emmap, scale_factors, radii)
    rp_after_shifted = compute_radial_profile(emmap_shifted_rp)
    freq = frequency_array(rp_from_emmap, 1.2156)
    plot_radial_profile(freq,[rp_from_emmap, rp_to_emmap, new_rp_from_emmap,rp_after_shifted])
    
    
    return emmap_shifted_rp
    
##### SAVE FUNCTIONS #####

def save_list_as_map(values_list, masked_indices, map_shape, map_path, apix):
    from locscale.include.emmer.ndimage.map_utils import save_as_mrc
    from locscale.utils.general import put_scaled_voxels_back_in_original_volume_including_padding
    value_map = put_scaled_voxels_back_in_original_volume_including_padding(values_list, masked_indices, map_shape)
    save_as_mrc(value_map, output_filename=map_path, apix=apix)

def write_out_final_volume_window_back_if_required(args, LocScaleVol, parsed_inputs_dict):
    from locscale.utils.general import pad_or_crop_volume
    from locscale.include.emmer.ndimage.map_utils import save_as_mrc
    from locscale.utils.plot_tools import make_locscale_report
    import mrcfile
    
    input_map = mrcfile.open(parsed_inputs_dict['emmap_path']).data
    
    wn = parsed_inputs_dict['wn']
    window_bleed_and_pad =parsed_inputs_dict['win_bleed_pad']
    apix = parsed_inputs_dict['apix']
        
    if window_bleed_and_pad:
        #map_shape = [(LocScaleVol.shape[0] - wn), (LocScaleVol.shape[1] - wn), (LocScaleVol.shape[2] - wn)]
        map_shape = input_map.shape
        LocScaleVol = pad_or_crop_volume(LocScaleVol, (map_shape))
    output_filename = args.outfile
    if args.dev_mode:
        output_filename = output_filename[:-4]+"_devmode.mrc"
    save_as_mrc(map_data=LocScaleVol, output_filename=output_filename, apix=apix, origin=0, verbose=True)
    
    if args.symmetry != "C1":
        resolution = parsed_inputs_dict['fsc_resolution']
        print("Imposing a symmetry condition of {}".format(args.symmetry))
        from locscale.include.symmetry_emda.symmetry_map import symmetrize_map_emda
        
        LocScaleVol_sym = symmetrize_map_emda(emmap_path=output_filename,pg=args.symmetry)

        output_filename = output_filename[:-4]+"_symmetrised.mrc"

        save_as_mrc(map_data=LocScaleVol_sym, output_filename=output_filename, apix=apix, origin=0, verbose=True)
    
    if args.output_report:
        make_locscale_report(args, parsed_inputs_dict, output_filename, window_bleed_and_pad)
    
    return LocScaleVol

##### MPI related functions #####

def split_sequence_evenly(seq, size):
    """
    >>> split_sequence_evenly(list(range(9)), 4)
    [[0, 1], [2, 3, 4], [5, 6], [7, 8]]
    >>> split_sequence_evenly(list(range(18)), 4)
    [[0, 1, 2, 3, 4], [5, 6, 7, 8], [9, 10, 11, 12, 13], [14, 15, 16, 17]]
    """
    from locscale.utils.math_tools import round_up_proper
    newseq = []
    splitsize = 1.0 / size * len(seq)
    for i in range(size):
        newseq.append(seq[round_up_proper(i * splitsize):round_up_proper((i+1) * splitsize)])
    return newseq

def merge_sequence_of_sequences(seq):
    """
    >>> merge_sequence_of_sequences([list(range(9)), list(range(3))])
    [0, 1, 2, 3, 4, 5, 6, 7, 8, 0, 1, 2]
    >>> merge_sequence_of_sequences([list(range(9)), [], list(range(3))])
    [0, 1, 2, 3, 4, 5, 6, 7, 8, 0, 1, 2]
    """
    newseq = [number for sequence in seq for number in sequence]

    return newseq

############################################################### 
################ CODE HELL ####################################

    #print("6) Printed locscale quality metrics")
    #quality_metrics, emmap_local_df, locscale_local_df = validation_metrics(args, parsed_input, locscale_path)
    #print("7) validation metrics done")
    #local_histogram_analysis_skew_kurt_emmap_fig = plot_regression(emmap_local_df, x_col="skew_emmap", y_col="kurtosis_emmap", 
    #                                                     title_text="Emmap local analysis: skew and kurtosis")
    #local_histogram_analysis_mean_var_emmap_fig = plot_linear_regression(emmap_local_df, x_col="mean_emmap", y_col="variance_emmap", 
    #                                                     title_text="Emmap local analysis: mean and variance")
    
    #local_histogram_analysis_skew_kurt_locscale_fig = plot_regression(locscale_local_df, x_col="skew_emmap", y_col="kurtosis_emmap", 
    #                                                     title_text="Locscale local analysis: skew and kurtosis")
    #local_histogram_analysis_mean_var_locscale_fig = plot_linear_regression(locscale_local_df, x_col="mean_emmap", y_col="variance_emmap", 
    #                                                     title_text="Locscale local analysis: mean and variance")
        
    
# def validation_metrics(args, parsed_inputs_dict, locscale_path):
#     print("Calculating map quality metrics...")
#     import pandas as pd
#     from locscale.utils.map_quality import map_quality_kurtosis, map_quality_pdb, local_histogram_analysis, measure_debye_pwlf
#     from locscale.include.emmer.ndimage.map_quality_tools import calculate_adjusted_surface_area
#     from locscale.preprocessing.headers import number_of_segments
#     ## Em-map based quality metrics
    
#     emmap_path = parsed_inputs_dict['emmap_path']
#     mask_path = parsed_inputs_dict['mask_path']
#     locscale_path = locscale_path
#     fsc_resolution = parsed_inputs_dict['fsc_resolution']
    
#     wilson_cutoff = parsed_inputs_dict['scale_factor_args']['wilson']
    
#     if args.model_coordinates is not None:
#         pdb_path = args.model_coordinates
#         calculate_pdb_based_metric = True
#     else:
#         calculate_pdb_based_metric = False
    
#     emmap_kurtosis = map_quality_kurtosis(emmap_path, mask_path)
#     locscale_kurtosis = map_quality_kurtosis(locscale_path, mask_path)
    
#     emmap_adjusted_surface_area = calculate_adjusted_surface_area(emmap_path, mask_path=mask_path, fsc_resolution=fsc_resolution)
#     locscale_adjusted_surface_area = calculate_adjusted_surface_area(locscale_path, mask_path=mask_path, fsc_resolution=fsc_resolution)
    
#     emmap_debye_slope = measure_debye_pwlf(emmap_path, wilson_cutoff, fsc_resolution, number_of_segments(fsc_resolution))
#     locscale_debye_slope = measure_debye_pwlf(locscale_path, wilson_cutoff, fsc_resolution, number_of_segments(fsc_resolution))
    
#     emmap_local_df, emmap_local_r_squared_skew_kurtosis, emmap_local_r_squared_mean_variance = local_histogram_analysis(emmap_path, mask_path, fsc_resolution)
#     locscale_local_df, locscale_local_r_squared_skew_kurtosis, locscale_local_r_squared_mean_variance = local_histogram_analysis(locscale_path, mask_path, 
#                                                                                                                                  fsc_resolution)
    
#     quality_metrics = {
#         'emmap_kurtosis': emmap_kurtosis,
#         'locscale_kurtosis': locscale_kurtosis,
#         'emmap_debye_slope':emmap_debye_slope,
#         'locscale_debye_slope':locscale_debye_slope,
#         'emmap_adjusted_surface_area':emmap_adjusted_surface_area,
#         'locscale_adjusted_surface_area':locscale_adjusted_surface_area,
#         'emmap_masked_local_r_squared_skew_kurtosis':emmap_local_r_squared_skew_kurtosis,
#         'emmap_masked_local_r_squared_mean_variance':emmap_local_r_squared_mean_variance,
#         'locscale_local_r_squared_skew_kurtosis':locscale_local_r_squared_skew_kurtosis,
#         'locscale_local_r_squared_mean_variance':locscale_local_r_squared_mean_variance}
    
#     if calculate_pdb_based_metric:
#         emmap_rscc_metric = map_quality_pdb(emmap_path, mask_path, pdb_path, test='rscc')
#         locscale_rscc_metric = map_quality_pdb(locscale_path, mask_path, pdb_path, test='rscc')
        
#         emmap_fsc_metric = map_quality_pdb(emmap_path, mask_path, pdb_path, test='fsc')
#         locscale_fsc_metric = map_quality_pdb(locscale_path, mask_path, pdb_path, test='fsc')
        
#         quality_metrics['emmap_rscc_metric'] = emmap_rscc_metric
#         quality_metrics['locscale_rscc_metric'] = locscale_rscc_metric
#         quality_metrics['emmap_fsc_metric'] = emmap_fsc_metric
#         quality_metrics['locscale_fsc_metric'] = locscale_fsc_metric
    
#     return quality_metrics, emmap_local_df, locscale_local_df

    
# def print_locscale_quality_metrics(parsed_inputs_dict, locscale_map):
#     from locscale.include.emmer.ndimage.profile_tools import compute_radial_profile, frequency_array, estimate_bfactor_through_pwlf
#     from scipy.stats import kurtosis
#     import mrcfile
#     from locscale.utils.map_quality import map_quality_kurtosis
    
#     ## Quality based on radial profile
#     rp_locscale_map = compute_radial_profile(locscale_map)
#     apix = parsed_inputs_dict['apix']
#     freq = frequency_array(rp_locscale_map, apix=apix)
#     wilson_cutoff = parsed_inputs_dict['scale_factor_args']['wilson']
#     fsc_cutoff = parsed_inputs_dict['fsc_resolution']
    
    
#     bfactor_locscale, _,(fit, z, slopes_locscale) = estimate_bfactor_through_pwlf(freq, rp_locscale_map, wilson_cutoff=wilson_cutoff, fsc_cutoff=fsc_cutoff)
    
#     breakpoints = (1/np.sqrt(z)).round(2)
#     debye_slope = slopes_locscale[1]
#     r_squared = fit.r_squared()
    
#     slopes_unsharp = [round(x,1) for x in parsed_inputs_dict['bfactor_info'][2]]
#     debye_slope_unsharp = slopes_unsharp[1]
#     bfactor_unsharp = parsed_inputs_dict['bfactor_info'][0]
    
#     ## Kurtosis metric
#     map_kurtosis = kurtosis(locscale_map.flatten())
#     emmap_unsharpened = mrcfile.open(parsed_inputs_dict['emmap_path']).data
#     unsharpened_kurtosis = map_quality_kurtosis(parsed_inputs_dict['emmap_path'], parsed_inputs_dict['mask_path'])
    
#     map_quality = {}
#     map_quality['kurtosis_unsharpened'] = round(unsharpened_kurtosis,2)
#     map_quality['kurtosis_locscale'] = round(map_kurtosis,2)
#     map_quality['debye_slope_unsharpened'] = debye_slope_unsharp
#     map_quality['debye_slope_locscale'] = round(debye_slope,1)
#     map_quality['bfactor_unsharp'] = bfactor_unsharp
#     map_quality['bfactor_locscale'] = bfactor_locscale
#     map_quality['pwlf_breakpoints'] = breakpoints
#     map_quality['pwlf_slopes'] = slopes_locscale.round(2)
#     map_quality['pwlf_r_sq'] = round(r_squared,2)
#     map_quality['Locscale-shape'] = locscale_map.shape
#     map_quality['Locscale-scale'] = [round(locscale_map.min(),2),round(locscale_map.max(),2)]
#     map_quality['Unsharpened_map-scale'] = [round(emmap_unsharpened.min(),2),round(emmap_unsharpened.max(),2)]
#     print("Map quality results: \n")
#     print(map_quality)
    
#     import matplotlib.pyplot as plt
#     fig, ax =plt.subplots(figsize=(16,16))
#     ax.axis('off')
    
#     text = []
#     for key in map_quality.keys():
#         text.append([key, map_quality[key]])
        
 
#     table= ax.table(cellText=text, loc="center", colLabels=["Parameter","Values"], cellLoc='center')
#     table.auto_set_font_size(False)
#     table.set_fontsize(16)
#     table.scale(1,4)
#     return fig



## MATH FUNCTIONS

# def round_up_proper(x):
#     epsilon = 1e-5  ## To round up in case of rounding to odd
#     return np.round(x+epsilon).astype(int)

# def round_up_to_even(x):
#     ceil_x = math.ceil(x)
#     if ceil_x % 2 == 0:   ## check if it's even, if not return one higher
#         return ceil_x
#     else:
#         return ceil_x+1

# def round_up_to_odd(x):
#     ceil_x = math.ceil(x)
#     if ceil_x % 2 == 0:   ## check if it's even, if so return one higher
#         return ceil_x+1
#     else:
#         return ceil_x

# def true_percent_probability(n):
#     x = np.random.uniform(low=0, high=100)
#     if x <= n:
#         return True
#     else:
#         return False

# def get_map_characteristics(parsed_inputs_dict, figsize=(12,12), font=12):
#     import matplotlib.pyplot as plt
    
#     fig, ax =plt.subplots(figsize=figsize)
    
#     ax.axis('off')
    
#     required_stats = {}
#     required_stats['UseTheoreticalProfiles'] = parsed_inputs_dict['use_theoretical']
#     required_stats['WindowSizePixel'] = parsed_inputs_dict['wn']
#     required_stats['apix'] = parsed_inputs_dict['apix']
#     required_stats['WindowBleedPad'] = parsed_inputs_dict['win_bleed_pad']
#     required_stats['EmmapShapeForLocScale'] = parsed_inputs_dict['emmap'].shape
#     required_stats['WilsonCutoff'] = round(parsed_inputs_dict['scale_factor_args']['wilson'],2)
#     required_stats['HighFreqCutoff'] = round(parsed_inputs_dict['scale_factor_args']['high_freq'],2)
#     required_stats['FSC'] = round(parsed_inputs_dict['scale_factor_args']['fsc_cutoff'],2)
#     required_stats['Bfactor'] = parsed_inputs_dict['bfactor_info'][0]
#     required_stats['Breakpoints'] = [round(x,1) for x in parsed_inputs_dict['bfactor_info'][1]]
#     required_stats['Slopes'] = [round(x,1) for x in parsed_inputs_dict['bfactor_info'][2]]
    
#     print("General statistics: \n")
#     print(required_stats)
    
#     text = []
#     for key in required_stats.keys():
#         text.append([key, required_stats[key]])
    
#     table= ax.table(cellText=text, loc="center", colLabels=["Parameter","Values"], cellLoc='center')
#     table.auto_set_font_size(False)
#     table.set_fontsize(font)
#     table.scale(1,4)
#     return fig

# def linear(x,a,b):
#     return a * x + b

# def general_quadratic(x,a,b,c):
#     return a * x**2 + b*x + c
    
# def r2(y_fit, y_data):
#     y_mean = y_data.mean()
#     residual_squares = (y_data-y_fit)**2
#     variance = (y_data-y_mean)**2
    
#     residual_sum_of_squares = residual_squares.sum()
#     sum_of_variance = variance.sum()
    
#     r_squared = 1 - residual_sum_of_squares/sum_of_variance
    
#     return r_squared

# ## FILE HANDLING FUNCTIONS
# def copy_file_to_folder(full_path_to_file, new_folder):
#     import shutil
#     import os
    
#     source = full_path_to_file
#     file_name = os.path.basename(source)
#     destination = os.path.join(new_folder, file_name)
#     shutil.copyfile(source, destination)
    
#     return destination

# def change_directory(args, folder_name):
#     import os    
#     from locscale.utils.file_tools import copy_file_to_folder
    
#     if folder_name == "processing_files":
#         current_directory = os.getcwd()
#         new_directory = os.path.join(current_directory, folder_name)
#     else:
#         new_directory = folder_name
    
#     if not os.path.isdir(new_directory):
#         os.mkdir(new_directory)
    
#     if args.verbose:
#         print("Copying files to {}\n".format(new_directory))

#     for arg in vars(args):
#         value = getattr(args, arg)
#         if isinstance(value, str):
#             if os.path.exists(value) and arg != "outfile" and arg != "output_processing_files":
#                 new_location=copy_file_to_folder(value, new_directory)
#                 setattr(args, arg, new_location)
    
#     if args.verbose:
#         print("Changing active directory to: {}\n".format(new_directory))
#     os.chdir(new_directory)
    
#     return args

# def generate_filename_from_halfmap_path(in_path):
#     ## find filename in the path    
#     filename = in_path.split("/")[-1]
    
#     ## find EMDB ID in filename
    
#     possible_emdb_id = [filename[x:x+4] for x in range(len(filename)-3) if filename[x:x+4].isnumeric()]
#     if len(possible_emdb_id) == 1:
#         emdb_id = possible_emdb_id[0]
#         newfilename = ["EMD_"+emdb_id+"_unfiltered.mrc"]
#     else:
#         newfilename = ["emdb_map_unfiltered.mrc"]
    
    
#     new_path = "/".join(in_path.split("/")[:-1]+newfilename)
    
#     return new_path
    
# def get_emmap_path_from_args(args):
#     from locscale.utils.file_tools import generate_filename_from_halfmap_path
#     from locscale.include.emmer.ndimage.map_tools import add_half_maps
    
#     if args.emmap_path is not None:    
#         emmap_path = args.emmap_path
#         shift_vector=shift_map_to_zero_origin(emmap_path)
#     elif args.halfmap_paths is not None:
#         print("Adding the two half maps provided to generate a full map \n")
#         halfmap_paths = args.halfmap_paths
#         assert len(halfmap_paths) == 2, "Please provide two half maps"

#         halfmap1_path = halfmap_paths[0]
#         halfmap2_path = halfmap_paths[1]
#         new_file_path = generate_filename_from_halfmap_path(halfmap1_path)
#         emmap_path = add_half_maps(halfmap1_path, halfmap2_path,new_file_path)
#         shift_vector=shift_map_to_zero_origin(halfmap1_path)
    
#     return emmap_path, shift_vector
        
#         ## TBC

# def is_input_path_valid(list_of_test_paths):
#     '''
#     Check if a list of paths are not None and if path points to an actual file

#     Parameters
#     ----------
#     list_of_test_paths : list
#         list of paths

#     Returns
#     -------
#     None.

#     '''
#     import os
    
#     for test_path in list_of_test_paths:
#         if test_path is None:
#             is_test_path_valid = False
#             return is_test_path_valid
#         if not os.path.exists(test_path):
#             is_test_path_valid = False
#             return is_test_path_valid
    
#     ## If all tests passed then return True
#     is_test_path_valid = True
#     return is_test_path_valid

# def check_user_input(args):
#     '''
#     Check user inputs for errors and conflicts

#     Parameters
#     ----------
#     args : TYPE
#         DESCRIPTION.

#     Returns
#     -------
#     None.

#     '''
#     return 
#     if args.dev_mode:
#         print("Warning: You are in Dev mode. Not checking user input! Results maybe unreliable")
#         return 
    
#     import mrcfile
#     print("ignore profiles default", args.ignore_profiles)
#     ## Check input files
#     emmap_absent = True
#     if args.emmap_path is not None:
#         if is_input_path_valid([args.emmap_path]):
#             emmap_absent = False
    
#     half_maps_absent = True
#     if args.half_map1 is not None and args.half_map2 is not None:
#         if is_input_path_valid([args.half_map1, args.half_map2]):
#             half_maps_absent = False
    
#     mask_absent = True
#     if args.mask is not None:
#         if is_input_path_valid([args.mask]):
#             mask_absent = False
    
#     model_map_absent = True
#     if args.model_map is not None:
#         if is_input_path_valid([args.model_map]):
#             model_map_absent = False
    
#     model_coordinates_absent = True
#     if args.model_coordinates is not None:
#         if is_input_path_valid([args.model_coordinates]):
#             model_coordinates_absent = False
    
#     ## Rename variables
#     emmap_present, half_maps_present = not(emmap_absent), not(half_maps_absent)
#     model_map_present, model_coordinates_present = not(model_map_absent), not(model_coordinates_absent)
#     ## Sanity checks
    
#     ## If emmap is absent or half maps are absent, raise Exceptions
    
#     if emmap_absent and half_maps_absent:
#         raise UserWarning("Please input either an unsharpened map or two half maps")
          
    
#     if model_coordinates_present and model_map_present:
#         raise UserWarning("Please provide either a model map or a model coordinates. Not both")
    
#     ## If neither model map or model coordinates are provided, then users cannot use --ignore_profiles and --skip_refine flags
#     if model_coordinates_absent and model_map_absent:
#         if args.ignore_profiles:
#             raise UserWarning("You have not provided a Model Map or Model Coordinates. Thus, pseudo-atomic model will be used for \
#                               local sharpening. Please do not raise the --ignore_profiles flag")
#         if args.skip_refine:
#             raise UserWarning("You have not provided a Model Map or Model Coordinates. Performing REFMAC refinement is essential for \
#                               succesful operation of the procedue. Please do not raise the --skip_refine flag")
        
#         if args.ref_resolution is None:
#             raise UserWarning("You have not provided a Model Map or Model Coordinates. To use REFMAC refinement, resolution target is necessary. \
#                               Please provide a target resolution using -res or --ref_resolution")
                            
    
#     if model_coordinates_present and model_map_absent:
#         if args.skip_refine:
#             print("You have asked to skip REFMAC refinement. Atomic bfactors from the input model will be used for simulating Model Map")
#         else:
#             if args.ref_resolution is None:
#                 raise UserWarning("You have provided Model Coordinates. By default, the model bfactors will be refined using REFMAC. \
#                                   For this, a target resolution is required. Please provide this resolution target using -res or --ref_resolution. \
#                                       Instead if you think model bfactors are accurate, then raise the --skip_refine flag to ignore bfactor refinement.")
            

   
    
#     ## Check for window size < 10 A
#     if args.window_size is not None:
#         window_size_pixels = int(args.window_size)
#         if window_size_pixels%2 > 0:
#             print("You have input an odd window size. For best performance, an even numbered window size is required. Adding 1 to the provided window size ")
#         if args.apix is not None:
#             apix = float(args.apix)
#         else:
#             if args.emmap_path is not None:
#                 apix = mrcfile.open(args.emmap_path).voxel_size.x
#             elif args.half_map1 is not None:
#                 apix = mrcfile.open(args.half_map1).voxel_size.x
        
#         window_size_ang = window_size_pixels * apix
        
        
#         if window_size_ang < 10:
#             print("Warning: Provided window size of {} is too small for pixel size of {}. \
#                   Default window size is generally 25 A. Think of increasing the window size".format(window_size_pixels, apix))
                  


#     if args.outfile is None:
#         print("You have not entered a filename for LocScale output. Using a standard output file name: loc_scale.mrc. \
#               Any file with the same name in the current directory will be overwritten")

#         outfile = [x for x in vars(args) if x=="outfile"]
        
#         setattr(args, outfile[0], "loc_scale.mrc")

# ## PLOT FUNCTIONS
# def plot_regression(data_input, x_col, y_col, x_label=None, y_label=None, title_text=None):
#     from matplotlib.offsetbox import AnchoredText
#     import matplotlib.pyplot as plt
#     from scipy.optimize import curve_fit
#     from locscale.utils.math_tools import general_quadratic, r2
    
    
#     f, ax = plt.subplots(1,1)

#     def get_sign(x, leading=False):
#         if x < 0:
#             return "-"
#         else:
#             if leading:
#                 return ""
#             else:
#                 return "+"
            
#     data_unsort = data_input.copy()
#     data=data_unsort.sort_values(by=x_col)
#     x_data = data[x_col]
#     y_data = data[y_col]
    
#     p_opt, p_cov = curve_fit(general_quadratic, x_data, y_data)
#     a,b,c = p_opt
    
#     y_fit = general_quadratic(x_data, *p_opt)
    
#     r_squared = r2(y_fit, y_data)
    
#     ax.plot(x_data, y_data,'bo')
#     ax.plot(x_data, y_fit, 'r-')
#     equation = "y = {} {} x$^2$ {} {} x {} {}".format(get_sign(a,True),round(abs(a),1),get_sign(b), round(abs(b),1),get_sign(c),round(abs(c),1))
#     legend_text = equation + "\n" + "R$^2$={}".format(round(r_squared,2))
#     anchored_text=AnchoredText(legend_text, loc=2)
#     ax.add_artist(anchored_text)
#     if x_label is not None:
#         ax.set_xlabel(x_label)
#     else:
#         ax.set_xlabel(x_col)
        
#     if y_label is not None:
#         ax.set_ylabel(y_label)
#     else:
#         ax.set_ylabel(y_col)
#     ax.set_title(title_text)
    
#     return f
    
# def plot_linear_regression(data_input, x_col, y_col, x_label=None, y_label=None, title_text=None):
#     from matplotlib.offsetbox import AnchoredText
#     import matplotlib.pyplot as plt
#     from scipy.optimize import curve_fit
#     from locscale.utils.math_tools import linear, r2
    
#     f, ax = plt.subplots(1,1)

#     def get_sign(x, leading=False):
#         if x < 0:
#             return "-"
#         else:
#             if leading:
#                 return ""
#             else:
#                 return "+"
            
#     data_unsort = data_input.copy()
#     data=data_unsort.sort_values(by=x_col)
#     x_data = data[x_col]
#     y_data = data[y_col]
    
#     p_opt, p_cov = curve_fit(linear, x_data, y_data)
#     a,b = p_opt
    
#     y_fit = linear(x_data, *p_opt)
    
#     r_squared = r2(y_fit, y_data)
    
#     ax.plot(x_data, y_data,'bo')
#     ax.plot(x_data, y_fit, 'r-')
#     equation = "y = {} {} x {} {} ".format(get_sign(a,True),round(abs(a),1),get_sign(b), round(abs(b),1))
#     legend_text = equation + "\n" + "R$^2$={}".format(round(r_squared,2))
#     anchored_text=AnchoredText(legend_text, loc=2)
#     ax.add_artist(anchored_text)
#     if x_label is not None:
#         ax.set_xlabel(x_label)
#     else:
#         ax.set_xlabel(x_col)
        
#     if y_label is not None:
#         ax.set_ylabel(y_label)
#     else:
#         ax.set_ylabel(y_col)
#     ax.set_title(title_text)  
    
#     return f

# def print_input_arguments(args):
#     import matplotlib.pyplot as plt
    
#     fig, ax =plt.subplots(figsize=(16,16))
    
#     ax.axis('off')
    
#     text = []
#     path_arguments = [x for x in vars(args) if x in ["emmap_path","half_map1","half_map2","model_map",
#                                                   "mask","model_coordinates","outfile"]]
#     for arg in vars(args):
#         val = getattr(args, arg)
#         if arg in path_arguments and val is not None:
#             full_path = val
#             filename = full_path.split("/")[-1]
#             text.append([arg, filename])
#         else:
#             text.append([arg, val])
    
    
#     table= ax.table(cellText=text, loc="center", colLabels=["Parameter","Values"], cellLoc='center')
#     table.auto_set_font_size(False)
#     table.set_fontsize(16)
#     table.scale(1,2)
#     return fig
   
# def plot_bfactor_distribution_standard(unsharpened_emmap_path, locscale_map_path, mask_path, fsc_resolution):
#     from locscale.include.emmer.ndimage.map_tools import get_bfactor_distribution
#     import matplotlib.pyplot as plt
#     import seaborn as sns
    
#     fig, ax =plt.subplots(figsize=(8,8))
    
#     unsharped_emmap_dist = get_bfactor_distribution(unsharpened_emmap_path, mask_path, fsc_resolution)
   
#     locscale_dist = get_bfactor_distribution(locscale_map_path, mask_path, fsc_resolution)
    
#     unsharpened_array = np.array([x[0] for x in unsharped_emmap_dist.values()])
#     locscale_array = np.array([x[0] for x in locscale_dist.values()])
   
#     sns.kdeplot(unsharpened_array)
   
#     sns.kdeplot(locscale_array)
    
#     plt.legend(["unsharpened map","Locscale map"])
#     return fig

# def plot_pickle_output(folder):
#     import pickle
#     import random
#     from locscale.include.emmer.ndimage.profile_tools import plot_radial_profile
#     import os
    
#     pickle_output = os.path.join(folder, "profiles_audit.pickle")
#     with open(pickle_output,"rb") as audit_file:
#         audit_scaling = pickle.load(audit_file)
    

#     random_positions = list(audit_scaling.keys())    
#     key = random.choice(random_positions)
    
#     freq = audit_scaling[key]['freq']
#     em_profile = audit_scaling[key]['em_profile']
#     ref_profile = audit_scaling[key]['input_ref_profile']
#     theoretical_profile = audit_scaling[key]['theoretical_amplitude']
#     scaled_theoretical = audit_scaling[key]['scaled_theoretical_amplitude']
#     deviated_profile = audit_scaling[key]['deviated_reference_profile']
#     exponential_fit = audit_scaling[key]['exponential_fit']
    
        
        
#     fig=plot_radial_profile(freq,[em_profile, ref_profile, theoretical_profile, scaled_theoretical, deviated_profile, exponential_fit],legends=['em_profile','ref_profile','th profile','scaled th profile','Deviated profile','exponential fit'])
    
#     return fig

# def make_locscale_report(args, parsed_input, locscale_path, window_bleed_and_pad, report_output_filename=None, statistic_output_filename=None):
#     from locscale.include.emmer.ndimage.profile_tools import plot_emmap_section
#     from locscale.include.emmer.ndimage.profile_tools import plot_radial_profile, compute_radial_profile, frequency_array 
#     from matplotlib.backends.backend_pdf import PdfPages
#     import os
#     import mrcfile
#     from locscale.include.emmer.ndimage.fsc_util import plot_fsc_maps
#     from locscale.utils.general import pad_or_crop_volume
    
#     ## Input-Output characteristics
#     locscale_map = mrcfile.open(locscale_path).data
    
    
#     processing_files_folder = parsed_input['processing_files_folder']
#     pdffile = os.path.join(processing_files_folder, args.report_filename+"_general.pdf")
#     pdf = PdfPages(pdffile)
    
#     print("Preparing LocScale report: \n {}".format(pdffile))
    
#     if window_bleed_and_pad:
#         from locscale.utils.general import pad_or_crop_volume
#         emmap = pad_or_crop_volume(parsed_input['emmap'], locscale_map.shape)
#         modmap = pad_or_crop_volume(parsed_input['modmap'], locscale_map.shape)
#     else:
#         emmap = parsed_input['emmap']
#         modmap = parsed_input['modmap']
  
#     rp_emmap = compute_radial_profile(emmap)
#     rp_modmap = compute_radial_profile(modmap)
#     rp_locscale = compute_radial_profile(locscale_map)
#     freq = frequency_array(rp_emmap, apix=parsed_input['apix'])
    
    
#     #1  Input Table
#     try:
#         input_table = print_input_arguments(args)
#         pdf.savefig(input_table)
#     except Exception as e:
#         print("Could not print input table ! \n")
#         print(e)
    
#     #2 Radial Profiles

#     try:
#         radial_profile_fig = plot_radial_profile(freq, [rp_emmap, rp_modmap, rp_locscale],legends=['input_emmap', 'model_map','locscale_map'])
#         pdf.savefig(radial_profile_fig)
#     except Exception as e:
#         print("Could not print radial profiles")
#         print(e)
        
#     #3 Sections
    
#     try:
#         emmap_section_fig = plot_emmap_section(parsed_input['emmap'], title="Input")
#         pdf.savefig(emmap_section_fig)
#     except Exception as e:
#         print("Could not print Emmap section")
#         print(e)
    
#     try:
#         locscale_section_fig = plot_emmap_section(locscale_map, title="LocScale Output")
#         pdf.savefig(locscale_section_fig)
#     except Exception as e:
#         print("Could not print Locscale map section")
#         print(e)
        
#     #4 FSC curves
    
        
#     try:
#         fsc_figure = plot_fsc_maps(emmap, locscale_map, apix=parsed_input['apix'], title="FSC curve unsharpened input and sharpened map", font=12)
#         pdf.savefig(fsc_figure)
#     except Exception as e:
#         print("Could not print fsc_figure")
#         print(e)
    
#     #5 Bfactor distributions
#     try:
#         bfactor_kde_fig = plot_bfactor_distribution_standard(unsharpened_emmap_path=parsed_input['emmap_path'],
#                                                  mask_path=parsed_input['mask_path'], locscale_map_path=locscale_path, fsc_resolution=parsed_input['fsc_resolution'])
#         pdf.savefig(bfactor_kde_fig)
#     except Exception as e:
#         print("Could not print bfactor_kde_fig")
#         print(e)
        
#     try:
#         stats_table = get_map_characteristics(parsed_input)
#         pdf.savefig(stats_table)
#     except Exception as e:
#         print("Could not print stats_table")
#         print(e)
  
#     try:      
#         if parsed_input['use_theoretical']:
#             pickle_output_sample_fig = plot_pickle_output(processing_files_folder)
#             pdf.savefig(pickle_output_sample_fig)
#     except Exception as e:
#         print("Could not print pickle_output_sample_fig")
#         print(e)
                
    
#     pdf.close()  
