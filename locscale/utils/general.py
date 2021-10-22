import numpy as np
import math


def round_up_to_even(x):
    ceil_x = math.ceil(x)
    if ceil_x % 2 == 0:   ## check if it's even, if not return one higher
        return ceil_x
    else:
        return ceil_x+1

def round_up_to_odd(x):
    ceil_x = math.ceil(x)
    if ceil_x % 2 == 0:   ## check if it's even, if so return one higher
        return ceil_x+1
    else:
        return ceil_x

def true_percent_probability(n):
    x = np.random.uniform(low=0, high=100)
    if x <= n:
        return True
    else:
        return False

def copy_file_to_folder(full_path_to_file, new_folder):
    import shutil
    
    source = full_path_to_file
    file_name = source.split("/")[-1]
    destination = "/".join(new_folder.split("/")+[file_name])
    shutil.copyfile(source, destination)
    
    return destination

def change_directory(args, folder_name="processed"):
    import os    
    from locscale.utils.general import copy_file_to_folder
    
    current_directory = os.getcwd()
    new_directory = "/".join(current_directory.split("/")+[folder_name])
    if not os.path.isdir(new_directory):
        os.mkdir(new_directory)
    
    for arg in vars(args):
        value = getattr(args, arg)
        if isinstance(value, str):
            if os.path.exists(value) and arg != "outfile":
                print("Copying {} to a new directory: \n{}".format(value, new_directory))
                new_location=copy_file_to_folder(value, new_directory)
                print("file saved at: {}".format(new_location))
                setattr(args, arg, new_location)
    
    os.chdir(new_directory)
    
    return args

def gather_statistics(parsed_inputs_dict):
    import matplotlib.pyplot as plt
    
    fig, ax =plt.subplots(figsize=(16,16))
    
    ax.axis('off')
    
    required_stats = {}
    required_stats['UseTheoreticalProfiles'] = parsed_inputs_dict['use_theoretical']
    required_stats['WindowSizePixel'] = parsed_inputs_dict['wn']
    required_stats['apix'] = parsed_inputs_dict['apix']
    required_stats['WindowBleedPad'] = parsed_inputs_dict['win_bleed_pad']
    required_stats['EmmapShapeForLocScale'] = parsed_inputs_dict['emmap'].shape
    required_stats['WilsonCutoff'] = round(parsed_inputs_dict['scale_factor_args']['wilson'],2)
    required_stats['HighFreqCutoff'] = round(parsed_inputs_dict['scale_factor_args']['high_freq'],2)
    required_stats['FSC'] = round(parsed_inputs_dict['scale_factor_args']['fsc_cutoff'],2)
    required_stats['Smooth'] = parsed_inputs_dict['scale_factor_args']['smooth']
    required_stats['Bfactor'] = parsed_inputs_dict['bfactor_info'][0]
    required_stats['Breakpoints'] = [round(x,1) for x in parsed_inputs_dict['bfactor_info'][1]]
    required_stats['Slopes'] = [round(x,1) for x in parsed_inputs_dict['bfactor_info'][2]]
    
    text = []
    for key in required_stats.keys():
        text.append([key, required_stats[key]])
        
 
    table= ax.table(cellText=text, loc="center", colLabels=["Parameter","Values"], cellLoc='center')
    table.auto_set_font_size(False)
    table.set_fontsize(16)
    table.scale(1,4)
    return fig

def print_locscale_quality_metrics(parsed_inputs_dict, locscale_map):
    from locscale.include.emmer.ndimage.profile_tools import compute_radial_profile, frequency_array, estimate_bfactor_through_pwlf
    from scipy.stats import kurtosis
    import mrcfile
    
    
    ## Quality based on radial profile
    rp_locscale_map = compute_radial_profile(locscale_map)
    apix = parsed_inputs_dict['apix']
    freq = frequency_array(rp_locscale_map, apix=apix)
    wilson_cutoff = parsed_inputs_dict['scale_factor_args']['wilson']
    fsc_cutoff = parsed_inputs_dict['fsc_resolution']
    
    
    bfactor_locscale, _,(fit, z, slopes_locscale) = estimate_bfactor_through_pwlf(freq, rp_locscale_map, wilson_cutoff=wilson_cutoff, fsc_cutoff=fsc_cutoff)
    
    breakpoints = (1/np.sqrt(z)).round(2)
    debye_slope = slopes_locscale[1]
    r_squared = fit.r_squared()
    
    slopes_unsharp = [round(x,1) for x in parsed_inputs_dict['bfactor_info'][2]]
    debye_slope_unsharp = slopes_unsharp[1]
    bfactor_unsharp = parsed_inputs_dict['bfactor_info'][0]
    
    ## Kurtosis metric
    map_kurtosis = kurtosis(locscale_map.flatten())
    emmap_unsharpened = mrcfile.open(parsed_inputs_dict['emmap_path']).data
    unsharpened_kurtosis = kurtosis(emmap_unsharpened.flatten())
    
    map_quality = {}
    map_quality['kurtosis_unsharpened'] = unsharpened_kurtosis
    map_quality['kurtosis_locscale'] = map_kurtosis
    map_quality['debye_slope_unsharpened'] = debye_slope_unsharp
    map_quality['debye_slope_locscale'] = debye_slope
    map_quality['bfactor_unsharp'] = bfactor_unsharp
    map_quality['bfactor_locscale'] = bfactor_locscale
    map_quality['pwlf_breakpoints'] = breakpoints
    map_quality['pwlf_slopes'] = slopes_locscale.round(2)
    map_quality['pwlf_r_sq'] = round(r_squared,2)
    map_quality['general:shape'] = locscale_map.shape
    map_quality['general:scale'] = [round(locscale_map.min(),2),round(locscale_map.max(),2)]
    
    import matplotlib.pyplot as plt
    fig, ax =plt.subplots(figsize=(16,16))
    ax.axis('off')
    
    text = []
    for key in map_quality.keys():
        text.append([key, map_quality[key]])
        
 
    table= ax.table(cellText=text, loc="center", colLabels=["Parameter","Values"], cellLoc='center')
    table.auto_set_font_size(False)
    table.set_fontsize(16)
    table.scale(1,4)
    return fig

def print_input_arguments(args):
    import matplotlib.pyplot as plt
    
    fig, ax =plt.subplots(figsize=(16,16))
    
    ax.axis('off')
    
    text = []
    path_arguments = [x for x in vars(args) if x in ["em_map","half_map1","half_map2","model_map",
                                                  "mask","model_coordinates","outfile"]]
    for arg in vars(args):
        val = getattr(args, arg)
        if arg in path_arguments and val is not None:
            full_path = val
            filename = full_path.split("/")[-1]
            text.append([arg, filename])
        else:
            text.append([arg, val])
    
    
    table= ax.table(cellText=text, loc="center", colLabels=["Parameter","Values"], cellLoc='center')
    table.auto_set_font_size(False)
    table.set_fontsize(16)
    table.scale(1,2)
    return fig
   

def make_locscale_report(args, parsed_input, locscale_map, window_bleed_and_pad):
    from locscale.include.emmer.ndimage.profile_tools import plot_emmap_section
    from locscale.include.emmer.ndimage.profile_tools import plot_radial_profile, compute_radial_profile, frequency_array 
    from matplotlib.backends.backend_pdf import PdfPages
    import os
    ## Input-Output characteristics
    
    cwd = os.getcwd()
    save_file_in_folder = "/".join(cwd.split("/")+["processing_files"])
    pdffile = "/".join(save_file_in_folder.split("/")+["locscale_report.pdf"])
    print("Preparing LocScale report: \n {}".format(pdffile))
    
    if window_bleed_and_pad:
        from locscale.utils.prepare_inputs import pad_or_crop_volume
        emmap = pad_or_crop_volume(parsed_input['emmap'], locscale_map.shape)
        modmap = pad_or_crop_volume(parsed_input['modmap'], locscale_map.shape)
    else:
        emmap = parsed_input['emmap']
        modmap = parsed_input['modmap']
        
    rp_emmap = compute_radial_profile(emmap)
    rp_modmap = compute_radial_profile(modmap)
    rp_locscale = compute_radial_profile(locscale_map)
    freq = frequency_array(rp_emmap, apix=parsed_input['apix'])

    radial_profile_fig = plot_radial_profile(freq, [rp_emmap, rp_modmap, rp_locscale],
                                             legends=['input_emmap', 'model_map','locscale_map'])
    
    emmap_section_fig = plot_emmap_section(parsed_input['emmap'], title="Input")
    
    locscale_section_fig = plot_emmap_section(locscale_map, title="LocScale Output")
    
    stats_table = gather_statistics(parsed_input)
    
    input_table = print_input_arguments(args)
    
    map_quality_table = print_locscale_quality_metrics(parsed_input, locscale_map)
    
    pdf = PdfPages(pdffile)
    pdf.savefig(input_table)
    pdf.savefig(radial_profile_fig)
    pdf.savefig(emmap_section_fig)
    pdf.savefig(locscale_section_fig)
    pdf.savefig(map_quality_table)
    pdf.savefig(stats_table)
    if parsed_input['use_theoretical']:
        pickle_output_sample_fig = plot_pickle_output(save_file_in_folder)
        pdf.savefig(pickle_output_sample_fig)
        
    
    pdf.close()
    
    
    

def plot_pickle_output(folder):
    import pickle
    import random
    from locscale.include.emmer.ndimage.profile_tools import plot_radial_profile
    
    pickle_output = "/".join(folder.split("/")+["profiles_audit.pickle"])
    with open(pickle_output,"rb") as audit_file:
        audit_scaling = pickle.load(audit_file)
    

    random_positions = list(audit_scaling.keys())    
    key = random.choice(random_positions)
    
    freq = audit_scaling[key]['freq']
    em_profile = audit_scaling[key]['em_profile']
    ref_profile = audit_scaling[key]['input_ref_profile']
    theoretical_profile = audit_scaling[key]['theoretical_amplitude']
    scaled_theoretical = audit_scaling[key]['scaled_theoretical_amplitude']
    merged_profile = audit_scaling[key]['scaled_reference_profile']
    
        
        
    fig=plot_radial_profile(freq,[em_profile, ref_profile, theoretical_profile, scaled_theoretical, merged_profile],legends=['em_profile','ref_profile','th profile','scaled th profile','merged'])
    
    return fig

def is_input_path_valid(list_of_test_paths):
    '''
    Check if a list of paths are not None and if path points to an actual file

    Parameters
    ----------
    list_of_test_paths : list
        list of paths

    Returns
    -------
    None.

    '''
    import os
    
    for test_path in list_of_test_paths:
        if test_path is None:
            is_test_path_valid = False
            return is_test_path_valid
        if not os.path.exists(test_path):
            is_test_path_valid = False
            return is_test_path_valid
    
    ## If all tests passed then return True
    is_test_path_valid = True
    return is_test_path_valid

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
    
    print("Current origin: ", current_origin)
    emmap_data = mrcfile.open(emmap_path).data
    
    output_file = emmap_path
    save_as_mrc(map_data=emmap_data, output_filename=emmap_path, apix=voxel_size, origin=0)
    
    shift_vector = target_origin - current_origin
    print("Shift vector: {} ".format(shift_vector.round(2)))
    return shift_vector


def round_up(x):
    return np.ceil(x).astype(int)

            
def get_emmap_path_from_args(args):
    from locscale.utils.prepare_inputs import generate_filename_from_halfmap_path
    from locscale.include.emmer.ndimage.map_tools import add_half_maps
    
    if args.em_map is not None:    
        emmap_path = args.em_map
        shift_vector=shift_map_to_zero_origin(args.em_map)
    elif args.half_map1 is not None and args.half_map2 is not None:
        print("Adding the two half maps provided to generate a full map \n")
        new_file_path = generate_filename_from_halfmap_path(args.half_map1)
        emmap_path = add_half_maps(args.half_map1, args.half_map2,new_file_path)
        shift_vector=shift_map_to_zero_origin(args.half_map1)
    
    return emmap_path, shift_vector
        
        ## TBC

def check_user_input(args):
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
    if args.god_mode:
        print("You are in God mode. Not checking user input!")
        return 
    
    import mrcfile
    print("ignore profiles default", args.ignore_profiles)
    ## Check input files
    emmap_absent = True
    if args.em_map is not None:
        if is_input_path_valid([args.em_map]):
            emmap_absent = False
    
    half_maps_absent = True
    if args.half_map1 is not None and args.half_map2 is not None:
        if is_input_path_valid([args.half_map1, args.half_map2]):
            half_maps_absent = False
    
    mask_absent = True
    if args.mask is not None:
        if is_input_path_valid([args.mask]):
            mask_absent = False
    
    model_map_absent = True
    if args.model_map is not None:
        if is_input_path_valid([args.model_map]):
            model_map_absent = False
    
    model_coordinates_absent = True
    if args.model_coordinates is not None:
        if is_input_path_valid([args.model_coordinates]):
            model_coordinates_absent = False
    
    ## Rename variables
    emmap_present, half_maps_present = not(emmap_absent), not(half_maps_absent)
    model_map_present, model_coordinates_present = not(model_map_absent), not(model_coordinates_absent)
    ## Sanity checks
    
    ## If emmap is absent or half maps are absent, raise Exceptions
    
    if emmap_absent and half_maps_absent:
        raise UserWarning("Please input either an unsharpened map or two half maps")
          
    
    if model_coordinates_present and model_map_present:
        raise UserWarning("Please provide either a model map or a model coordinates. Not both")
    
    ## If neither model map or model coordinates are provided, then users cannot use --ignore_profiles and --skip_refine flags
    if model_coordinates_absent and model_map_absent:
        if args.ignore_profiles:
            raise UserWarning("You have not provided a Model Map or Model Coordinates. Thus, pseudo-atomic model will be used for \
                              local sharpening. Please do not raise the --ignore_profiles flag")
        if args.skip_refine:
            raise UserWarning("You have not provided a Model Map or Model Coordinates. Performing REFMAC refinement is essential for \
                              succesful operation of the procedue. Please do not raise the --skip_refine flag")
        
        if args.ref_resolution is None:
            raise UserWarning("You have not provided a Model Map or Model Coordinates. To use REFMAC refinement, resolution target is necessary. \
                              Please provide a target resolution using -res or --ref_resolution")
                            
    
    if model_coordinates_present and model_map_absent:
        if args.skip_refine:
            print("You have asked to skip REFMAC refinement. Atomic bfactors from the input model will be used for simulating Model Map")
        else:
            if args.ref_resolution is None:
                raise UserWarning("You have provided Model Coordinates. By default, the model bfactors will be refined using REFMAC. \
                                  For this, a target resolution is required. Please provide this resolution target using -res or --ref_resolution. \
                                      Instead if you think model bfactors are accurate, then raise the --skip_refine flag to ignore bfactor refinement.")
            

   
    
    ## Check for window size < 10 A
    if args.window_size is not None:
        window_size_pixels = int(args.window_size)
        if args.apix is not None:
            apix = float(args.apix)
        else:
            if args.em_map is not None:
                apix = mrcfile.open(args.em_map).voxel_size.x
            elif args.half_map1 is not None:
                apix = mrcfile.open(args.half_map1).voxel_size.x
        
        window_size_ang = window_size_pixels * apix
        
        if window_size_ang < 10:
            print("Warning: Provided window size of {} is too small for pixel size of {}. \
                  Default window size is generally 25 A. Think of increasing the window size".format(window_size_pixels, apix))
                  


    if args.outfile is None:
        print("You have not entered a filename for LocScale output. Using a standard output file name: loc_scale.mrc. \
              Any file with the same name in the current directory will be overwritten")

        outfile = [x for x in vars(args) if x=="outfile"]
        
        setattr(args, outfile[0], "loc_scale.mrc")
    

        
    
    
