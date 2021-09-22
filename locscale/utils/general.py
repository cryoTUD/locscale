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
    import pandas as pd
    import matplotlib.pyplot as plt
    
    fig, ax =plt.subplots(figsize=(12,4))
    required_stats = {}
    required_stats['use_theoretical'] = parsed_inputs_dict['use_theoretical']
    required_stats['wilson'] = parsed_inputs_dict['scale_factor_args']['wilson']
    required_stats['high_freq'] = parsed_inputs_dict['scale_factor_args']['high_freq']
    required_stats['fsc_cutoff'] = parsed_inputs_dict['scale_factor_args']['fsc_cutoff']
    required_stats['smooth'] = parsed_inputs_dict['scale_factor_args']['smooth']
    required_stats['bfactor'] = parsed_inputs_dict['bfactor_info'][0]
    required_stats['breakpoints'] = parsed_inputs_dict['bfactor_info'][1]
    required_stats['slopes'] = parsed_inputs_dict['bfactor_info'][2]
    
    df = pd.DataFrame(data=required_stats.values(), index=required_stats.keys())
    ax.table(cellText=df.values, colLabels=df.index, loc="center")
    
    return fig
    

def make_locscale_report(parsed_input, locscale_map):
    from locscale.include.emmer.ndimage.profile_tools import plot_emmap_section
    from locscale.include.emmer.ndimage.profile_tools import plot_radial_profile, compute_radial_profile, frequency_array 
    import pandas as pd
    from matplotlib.backends.backend_pdf import PdfPages
    import os
    ## Input-Output characteristics
    
    cwd = os.getcwd()
    pdffile = "/".join(cwd.split("/")+["processing_files","locscale_report.pdf"])
    print("Preparing LocScale report: \n {}".format(pdffile))
    
    rp_emmap = compute_radial_profile(parsed_input['emmap'])
    rp_modmap = compute_radial_profile(parsed_input['modmap'])
    rp_locscale = compute_radial_profile(locscale_map)
    
    freq = frequency_array(rp_emmap, apix=parsed_input['apix'])
    
    radial_profile_fig = plot_radial_profile(freq, [rp_emmap, rp_modmap, rp_locscale],
                                             legends=['input_emmap', 'model_map','locscale_map'])
    
    emmap_section_fig = plot_emmap_section(parsed_input['emmap'])
    locscale_section_fig = plot_emmap_section(locscale_map)
    
    stats_table = gather_statistics(parsed_input)
    
    
    
    pdf = PdfPages(pdffile)
    pdf.savefig(radial_profile_fig)
    pdf.savefig(emmap_section_fig)
    pdf.savefig(locscale_section_fig)
    pdf.savefig(stats_table)
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
        