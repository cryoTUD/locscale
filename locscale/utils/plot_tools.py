
## PLOT FUNCTIONS
def plot_regression(data_input, x_col, y_col, x_label=None, y_label=None, title_text=None):
    from matplotlib.offsetbox import AnchoredText
    import matplotlib.pyplot as plt
    from scipy.optimize import curve_fit
    from locscale.utils.math_tools import general_quadratic, r2
    
    
    f, ax = plt.subplots(1,1)

    def get_sign(x, leading=False):
        if x < 0:
            return "-"
        else:
            if leading:
                return ""
            else:
                return "+"
            
    data_unsort = data_input.copy()
    data=data_unsort.sort_values(by=x_col)
    x_data = data[x_col]
    y_data = data[y_col]
    
    p_opt, p_cov = curve_fit(general_quadratic, x_data, y_data)
    a,b,c = p_opt
    
    y_fit = general_quadratic(x_data, *p_opt)
    
    r_squared = r2(y_fit, y_data)
    
    ax.plot(x_data, y_data,'bo')
    ax.plot(x_data, y_fit, 'r-')
    equation = "y = {} {} x$^2$ {} {} x {} {}".format(get_sign(a,True),round(abs(a),1),get_sign(b), round(abs(b),1),get_sign(c),round(abs(c),1))
    legend_text = equation + "\n" + "R$^2$={}".format(round(r_squared,2))
    anchored_text=AnchoredText(legend_text, loc=2)
    ax.add_artist(anchored_text)
    if x_label is not None:
        ax.set_xlabel(x_label)
    else:
        ax.set_xlabel(x_col)
        
    if y_label is not None:
        ax.set_ylabel(y_label)
    else:
        ax.set_ylabel(y_col)
    ax.set_title(title_text)
    
    return f
    
def plot_linear_regression(data_input, x_col, y_col, x_label=None, y_label=None, title_text=None):
    from matplotlib.offsetbox import AnchoredText
    import matplotlib.pyplot as plt
    from scipy.optimize import curve_fit
    from locscale.utils.math_tools import linear, r2
    
    f, ax = plt.subplots(1,1)

    def get_sign(x, leading=False):
        if x < 0:
            return "-"
        else:
            if leading:
                return ""
            else:
                return "+"
            
    data_unsort = data_input.copy()
    data=data_unsort.sort_values(by=x_col)
    x_data = data[x_col]
    y_data = data[y_col]
    
    p_opt, p_cov = curve_fit(linear, x_data, y_data)
    a,b = p_opt
    
    y_fit = linear(x_data, *p_opt)
    
    r_squared = r2(y_fit, y_data)
    
    ax.plot(x_data, y_data,'bo')
    ax.plot(x_data, y_fit, 'r-')
    equation = "y = {} {} x {} {} ".format(get_sign(a,True),round(abs(a),1),get_sign(b), round(abs(b),1))
    legend_text = equation + "\n" + "R$^2$={}".format(round(r_squared,2))
    anchored_text=AnchoredText(legend_text, loc=2)
    ax.add_artist(anchored_text)
    if x_label is not None:
        ax.set_xlabel(x_label)
    else:
        ax.set_xlabel(x_col)
        
    if y_label is not None:
        ax.set_ylabel(y_label)
    else:
        ax.set_ylabel(y_col)
    ax.set_title(title_text)  
    
    return f


class tab_print():
    def __init__(self, tab_size):
        self.tab_size = tab_size
        
    def tprint(self, string):
        text = "\t"*self.tab_size + string 
        print(text)
    


def print_input_arguments(args):
    import matplotlib.pyplot as plt
    
    fig, ax =plt.subplots(figsize=(16,16))
    
    ax.axis('off')
    
    text = []
    path_arguments = [x for x in vars(args) if x in ["emmap_path","half_map1","half_map2","model_map",
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
   
def plot_bfactor_distribution_standard(unsharpened_emmap_path, locscale_map_path, mask_path, fsc_resolution):
    from locscale.include.emmer.ndimage.map_tools import get_bfactor_distribution
    import matplotlib.pyplot as plt
    import seaborn as sns
    
    fig, ax =plt.subplots(figsize=(8,8))
    
    unsharped_emmap_dist = get_bfactor_distribution(unsharpened_emmap_path, mask_path, fsc_resolution)
   
    locscale_dist = get_bfactor_distribution(locscale_map_path, mask_path, fsc_resolution)
    
    unsharpened_array = np.array([x[0] for x in unsharped_emmap_dist.values()])
    locscale_array = np.array([x[0] for x in locscale_dist.values()])
   
    sns.kdeplot(unsharpened_array)
   
    sns.kdeplot(locscale_array)
    
    plt.legend(["unsharpened map","Locscale map"])
    return fig

def plot_pickle_output(folder):
    import pickle
    import random
    from locscale.include.emmer.ndimage.profile_tools import plot_radial_profile
    import os
    
    pickle_output = os.path.join(folder, "profiles_audit.pickle")
    with open(pickle_output,"rb") as audit_file:
        audit_scaling = pickle.load(audit_file)
    

    random_positions = list(audit_scaling.keys())    
    key = random.choice(random_positions)
    
    freq = audit_scaling[key]['freq']
    em_profile = audit_scaling[key]['em_profile']
    ref_profile = audit_scaling[key]['input_ref_profile']
    theoretical_profile = audit_scaling[key]['theoretical_amplitude']
    scaled_theoretical = audit_scaling[key]['scaled_theoretical_amplitude']
    deviated_profile = audit_scaling[key]['deviated_reference_profile']
    exponential_fit = audit_scaling[key]['exponential_fit']
    
        
        
    fig=plot_radial_profile(freq,[em_profile, ref_profile, theoretical_profile, scaled_theoretical, deviated_profile, exponential_fit],legends=['em_profile','ref_profile','th profile','scaled th profile','Deviated profile','exponential fit'])
    
    return fig

def make_locscale_report(args, parsed_input, locscale_path, window_bleed_and_pad, report_output_filename=None, statistic_output_filename=None):
    from locscale.include.emmer.ndimage.plots import plot_emmap_section, plot_radial_profile
    from locscale.utils.plot_tools import plot_pickle_output
    from locscale.include.emmer.ndimage.profile_tools import compute_radial_profile, frequency_array 
    from matplotlib.backends.backend_pdf import PdfPages
    import os
    import mrcfile
    from locscale.include.emmer.ndimage.fsc_util import plot_fsc_maps
    from locscale.utils.general import pad_or_crop_volume
    
    ## Input-Output characteristics
    locscale_map = mrcfile.open(locscale_path).data
    
    
    processing_files_folder = parsed_input['processing_files_folder']
    pdffile = os.path.join(processing_files_folder, args.report_filename+"_general.pdf")
    pdf = PdfPages(pdffile)
    
    print("Preparing LocScale report: \n {}".format(pdffile))
    
    if window_bleed_and_pad:
        from locscale.utils.general import pad_or_crop_volume
        emmap = pad_or_crop_volume(parsed_input['emmap'], locscale_map.shape)
        modmap = pad_or_crop_volume(parsed_input['modmap'], locscale_map.shape)
    else:
        emmap = parsed_input['emmap']
        modmap = parsed_input['modmap']
  
    rp_emmap = compute_radial_profile(emmap)
    rp_modmap = compute_radial_profile(modmap)
    rp_locscale = compute_radial_profile(locscale_map)
    freq = frequency_array(rp_emmap, apix=parsed_input['apix'])
    
    
    #1  Input Table
    try:
        input_table = print_input_arguments(args)
        pdf.savefig(input_table)
    except Exception as e:
        print("Could not print input table ! \n")
        print(e)
    
    #2 Radial Profiles

    try:
        radial_profile_fig = plot_radial_profile(freq, [rp_emmap, rp_modmap, rp_locscale],legends=['input_emmap', 'model_map','locscale_map'])
        pdf.savefig(radial_profile_fig)
    except Exception as e:
        print("Could not print radial profiles")
        print(e)
        
    #3 Sections
    
    try:
        emmap_section_fig = plot_emmap_section(parsed_input['emmap'], title="Input")
        pdf.savefig(emmap_section_fig)
    except Exception as e:
        print("Could not print Emmap section")
        print(e)
    
    try:
        locscale_section_fig = plot_emmap_section(locscale_map, title="LocScale Output")
        pdf.savefig(locscale_section_fig)
    except Exception as e:
        print("Could not print Locscale map section")
        print(e)
        
    #4 FSC curves
    
        
    try:
        fsc_figure = plot_fsc_maps(emmap, locscale_map, apix=parsed_input['apix'], title="FSC curve unsharpened input and sharpened map", font=12)
        pdf.savefig(fsc_figure)
    except Exception as e:
        print("Could not print fsc_figure")
        print(e)
    
    #5 Bfactor distributions
    try:
        bfactor_kde_fig = plot_bfactor_distribution_standard(unsharpened_emmap_path=parsed_input['emmap_path'],
                                                 mask_path=parsed_input['mask_path'], locscale_map_path=locscale_path, fsc_resolution=parsed_input['fsc_resolution'])
        pdf.savefig(bfactor_kde_fig)
    except Exception as e:
        print("Could not print bfactor_kde_fig")
        print(e)
        
    try:
        stats_table = get_map_characteristics(parsed_input)
        pdf.savefig(stats_table)
    except Exception as e:
        print("Could not print stats_table")
        print(e)
  
    try:      
        if parsed_input['use_theoretical']:
            pickle_output_sample_fig = plot_pickle_output(processing_files_folder)
            pdf.savefig(pickle_output_sample_fig)
    except Exception as e:
        print("Could not print pickle_output_sample_fig")
        print(e)
                
    
    pdf.close()  
