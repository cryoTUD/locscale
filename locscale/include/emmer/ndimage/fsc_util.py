import numpy as np
#from fsc_attributes_calc import *
#from spline_fit_fsc import *

def calculate_edge(wn):
    '''
    calculate edge for smoothing
    '''
    edge = min(int(wn/2)-3,6)
    return edge

def make_soft_edged_window(wn_shape,edge=5):
    z,y,x = wn_shape
#     z,y,x = (8,8,8)
#     edge = 2
    radius = int(round(max(z,y,x) / 2.0))
    rad_z = np.arange(np.floor(z/2.0)*-1, 
                      np.ceil(z/2.0))
    rad_y = np.arange(np.floor(y/2.0)*-1, 
                      np.ceil(y/2.0))
    rad_x = np.arange(np.floor(x/2.0)*-1, 
                      np.ceil(x/2.0))
    rad_x = rad_x**2
    rad_y = rad_y**2
    rad_z = rad_z**2
    dist = np.sqrt(rad_z[:,None,None]+rad_y[:,None] + rad_x)
    #for tanh smoothing get values from center as >=2 to -1 at fixed_radius
    fixed_radius = radius - edge
    fixed_radius = max(3,fixed_radius)
    dist = fixed_radius-dist
    dist[:] = (np.tanh(dist)+1)/2.
    dist[dist<0.] = 0.
    dist[dist>1.] = 1.
    return dist

def compare_tuple(tuple1,tuple2):
    for val1, val2 in zip(tuple1, tuple2):
        if type(val2) is float:
            if round(val1,2) != round(val2,2):
                return False
        else:
            if val1 != val2:
                return False
    return True

def calculate_shell_correlation(shell1,shell2):
    '''
    Calculate FSC in a resolution shell
    '''
    cov_ps1_ps2 = shell1*np.conjugate(shell2)
    sig_ps1 = shell1*np.conjugate(shell1)
    sig_ps2 = shell2*np.conjugate(shell2)
    cov_ps1_ps2 = np.sum(np.real(cov_ps1_ps2))
    var_ps1 = np.sum(np.real(sig_ps1))
    var_ps2 = np.sum(np.real(sig_ps2))
    #skip shells with no variance
    if np.round(var_ps1,15) == 0.0 or np.round(var_ps2,15) == 0.0: 
        fsc = 0.0
    else: fsc = cov_ps1_ps2/(np.sqrt(var_ps1*var_ps2))
    return fsc
    
def calculate_fsc(ps1,ps2,radii,map_shape):
    '''
    Calculate FSC curve given two FTmaps
    '''
    list_fsc = []
    list_radii = []
    list_nsf = []
    for r in np.unique(radii)[0:map_shape[0]//2]:
        idx = radii == r
        fsc = calculate_shell_correlation(ps1[idx],ps2[idx])
        list_fsc.append(fsc)
        list_radii.append(float(r)/(map_shape[0]))
        num_nonzero_avg = \
                min(np.count_nonzero(ps1[idx]),np.count_nonzero(ps2[idx]))
        list_nsf.append(num_nonzero_avg)
    #return resolution @ FSC 0.5
    if list_fsc[0] == -1.:
        list_fsc[0] = 1.
    list_fsc[0] = max(0.,list_fsc[0])
    #sorted tuples
    listfreq, listfsc, listnsf = zip(*sorted(zip(list_radii, list_fsc, list_nsf))) 
    return listfreq,listfsc,listnsf

def calculate_phase_correlation(ps1,ps2,radii,map_shape):
    '''
    Calculate FSC curve given two FTmaps
    '''
    list_fsc = []
    list_radii = []
    list_nsf = []
    ps1_phase = np.angle(ps1)
    ps2_phase = np.angle(ps2)
    for r in np.unique(radii)[0:map_shape[0]//2]:
        idx = radii == r
        fsc = calculate_shell_correlation(ps1_phase[idx],ps2_phase[idx])
        list_fsc.append(fsc)
        list_radii.append(float(r)/(map_shape[0]))
        num_nonzero_avg = \
                min(np.count_nonzero(ps1_phase[idx]),np.count_nonzero(ps2_phase[idx]))
        list_nsf.append(num_nonzero_avg)
    #return resolution @ FSC 0.5
    if list_fsc[0] == -1.:
        list_fsc[0] = 1.
    list_fsc[0] = max(0.,list_fsc[0])
    #sorted tuples
    listfreq, listfsc, listnsf = zip(*sorted(zip(list_radii, list_fsc, list_nsf))) 
    return listfreq,listfsc,listnsf

def calculate_amplitude_correlation(ps1,ps2,radii,map_shape):
    '''
    Calculate FSC curve given two FTmaps
    '''
    list_fsc = []
    list_radii = []
    list_nsf = []
    ps1_abs = np.abs(ps1)
    ps2_abs = np.abs(ps2)
    for r in np.unique(radii)[0:map_shape[0]//2]:
        idx = radii == r
        fsc = calculate_shell_correlation(ps1_abs[idx],ps2_abs[idx])
        list_fsc.append(fsc)
        list_radii.append(float(r)/(map_shape[0]))
        num_nonzero_avg = \
                min(np.count_nonzero(ps1_abs[idx]),np.count_nonzero(ps2_abs[idx]))
        list_nsf.append(num_nonzero_avg)
    #return resolution @ FSC 0.5
    if list_fsc[0] == -1.:
        list_fsc[0] = 1.
    list_fsc[0] = max(0.,list_fsc[0])
    #sorted tuples
    listfreq, listfsc, listnsf = zip(*sorted(zip(list_radii, list_fsc, list_nsf))) 
    return listfreq,listfsc,listnsf

def plot_fscs_old(dict_points,outfile,xlabel=None,ylabel=None,map_apix=1.0,
              xlim=None,ylim=None,line=True,marker=True,lstyle=True,
              maxRes=1.5,minRes=20.0):
    '''
    Plot fscs given multiple lists of freq and fscs in a dictionary
    Area between minRes and maxRes will be shaded
    '''
    try:
        import matplotlib.pyplot as plt
    except RuntimeError:
        plt = None
    try: plt.style.use('ggplot')
    except AttributeError: pass
    ymaxm = xmaxm = -100.0
    for k in dict_points:
        if max(dict_points[k][1]) > ymaxm: ymaxm = max(dict_points[k][1])
        if max(dict_points[k][0]) > xmaxm: xmaxm = max(dict_points[k][0]) 
        colormap = plt.cm.brg#Set1,Spectral#YlOrRd,Spectral,BuGn,Set1,Accent,spring
    if len(dict_points) < 4: colormap = plt.cm.gist_earth
    plt.gca().set_prop_cycle('color',[colormap(i) for i in np.linspace(0, 1, len(dict_points)+1)])
    if ylim is not None: plt.gca().set_ylim(ylim)
    plt.rcParams.update({'font.size': 18})
    plt.rcParams.update({'legend.fontsize': 14})
        
    if not xlabel is None: plt.xlabel(xlabel, fontsize=15)
    if not ylabel is None: plt.ylabel(ylabel,fontsize=15)
    list_styles = []
    for i in range(0,len(dict_points),4):
        list_styles.extend(['-',':','-.','--'])
    list_markers = ['o', '*','>','D','s','p','<','v',':','h','x','+',',','.','_','2','d','^', 'H']
    while len(list_markers) < len(dict_points):
        list_markers.extend(list_markers)
    i = 0
    
    for k in dict_points:
        if line and marker: plt.plot(dict_points[k][0],dict_points[k][1],
                                     linewidth=2.0,label=k,
                                     linestyle=list_styles[i],
                                     marker=list_markers[i])
        elif line and lstyle: plt.plot(dict_points[k][0],dict_points[k][1],
                                       linewidth=2.0,label=k,
                                       linestyle=list_styles[i])
        elif line: plt.plot(dict_points[k][0],dict_points[k][1],
                            linewidth=1.0,label=k,color='g')
        elif marker: plt.plot(dict_points[k][0],dict_points[k][1],
                              label=k,marker=list_markers[i])
        i += 1
        if i == 1:
            x_array = dict_points[k][0]
            
    plt.axhline(y=0.5, color='black',linestyle='--')
    plt.axvspan(map_apix/minRes,map_apix/maxRes,alpha=0.5,color='grey')
    # Set the ticks and labels...
    locs,labs = plt.xticks()
    step = (max(locs)-min(locs))/10.
    locs = np.arange(min(locs),max(locs)+step,step)
    labels = np.round(map_apix/locs[1:],1)
    plt.xticks(locs[1:], labels,rotation='vertical')
    plt.savefig(outfile)
    plt.close()

def calculate_fsc_maps(input_map_1, input_map_2):
    '''
    Wrapper to calculate FSC curve from the above functions

    Parameters
    ----------
    input_map_1 : TYPE
        DESCRIPTION.
    input_map_2 : TYPE
        DESCRIPTION.

    Returns
    -------
    None.

    '''
    
    from locscale.include.emmer.ndimage.profile_tools import compute_radial_profile, frequency_array
    from locscale.include.emmer.ndimage.map_utils import parse_input
    
    emmap_1 = parse_input(input_map_1)
    emmap_2 = parse_input(input_map_2)
    
    fft_1 = np.fft.rfftn(emmap_1)
    fft_2 = np.fft.rfftn(emmap_2)
    
    _, radii = compute_radial_profile(emmap_1, return_indices=True)
    
    map_shape = emmap_1.shape
    
    _, fsc, _ = calculate_fsc(fft_1, fft_2, radii, map_shape)
    fsc = np.array(fsc)
    return np.array(fsc)

def calculate_phase_correlation_maps(input_map_1, input_map_2):
    '''
    Wrapper to calculate FSC curve from the above functions

    Parameters
    ----------
    input_map_1 : TYPE
        DESCRIPTION.
    input_map_2 : TYPE
        DESCRIPTION.

    Returns
    -------
    None.

    '''
    
    from locscale.include.emmer.ndimage.profile_tools import compute_radial_profile, frequency_array
    from locscale.include.emmer.ndimage.map_utils import parse_input
    
    emmap_1 = parse_input(input_map_1)
    emmap_2 = parse_input(input_map_2)
    
    fft_1 = np.fft.rfftn(emmap_1)
    fft_2 = np.fft.rfftn(emmap_2)
    
    _, radii = compute_radial_profile(emmap_1, return_indices=True)
    
    map_shape = emmap_1.shape
    
    _, fsc, _ = calculate_phase_correlation(fft_1, fft_2, radii, map_shape)
    fsc = np.array(fsc)
    return np.array(fsc)

def calculate_amplitude_correlation_maps(input_map_1, input_map_2):
    '''
    Wrapper to calculate FSC curve from the above functions

    Parameters
    ----------
    input_map_1 : TYPE
        DESCRIPTION.
    input_map_2 : TYPE
        DESCRIPTION.

    Returns
    -------
    None.

    '''
    
    from locscale.include.emmer.ndimage.profile_tools import compute_radial_profile, frequency_array
    from locscale.include.emmer.ndimage.map_utils import parse_input
    
    emmap_1 = parse_input(input_map_1)
    emmap_2 = parse_input(input_map_2)
    
    fft_1 = np.fft.rfftn(emmap_1)
    fft_2 = np.fft.rfftn(emmap_2)
    
    _, radii = compute_radial_profile(emmap_1, return_indices=True)
    
    map_shape = emmap_1.shape
    
    _, fsc, _ = calculate_amplitude_correlation(fft_1, fft_2, radii, map_shape)
    fsc = np.array(fsc)
    return np.array(fsc)


def plot_fscs(freq, list_of_fsc, colors=['r','g','b','k','y','m'], legends=None, font=12,showlegend=True, showPoints=True):
    import matplotlib.pyplot as plt
    
    if len(list_of_fsc) > 6:
        print("Enter maximum of 6 profiles only if you want to see colors")
        
    i = 0
    
    if showPoints:
        colors = [x+".-" for x in colors]
        
    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    ax1.grid(True)
    ax2 = ax1.twiny()
    for fsc in list_of_fsc:
        ax1.plot(freq,fsc, colors[i])
        i = i+1
            
    ax2.set_xticks(ax1.get_xticks())
    ax2.set_xbound(ax1.get_xbound())
    ax2.set_xticklabels([round(1/x,1) for x in ax1.get_xticks()])
    if legends is None:
        legends = ["FSC_"+str(i) for i in range(len(list_of_fsc))]
        
        
    ax1.legend(legends,fontsize=font)
    ax1.set_xlabel(r'$1/d [\AA^{-1}]$',fontsize=font)
    ax1.set_ylabel('FSC',fontsize=font)
    ax2.set_xlabel('$d [\AA]$',fontsize=font)

def plot_multiple_fsc(list_of_halfmap_tuple, common_mask_path, softening_parameter=5, apix=None, legend=None, title=None, font=16, colors=None):
    import matplotlib.pyplot as plt
    from matplotlib.pyplot import cm
    from scipy.interpolate import interp1d
    from locscale.include.emmer.ndimage.map_utils import parse_input
    from locscale.include.emmer.ndimage.profile_tools import frequency_array
    from locscale.include.emmer.ndimage.filter import get_cosine_mask
    import mrcfile
    from tqdm import tqdm
    
    mask = parse_input(common_mask_path)
    if apix is not None:
        apix = apix
    else:
        apix = mrcfile.open(common_mask_path).voxel_size.tolist()[0]
    
    softmask = get_cosine_mask(mask, length_cosine_mask_1d=softening_parameter)
    
    fsc_curves = {}
    
    for i,halfmap_tuple in enumerate(tqdm(list_of_halfmap_tuple,desc="Computing FSCs")):
        
        halfmap1 = parse_input(halfmap_tuple[0])
        halfmap2 = parse_input(halfmap_tuple[1])
        
        masked_halfmap1 = softmask * halfmap1
        masked_halfmap2 = softmask * halfmap2
        
        fsc_curve = calculate_fsc_maps(input_map_1=masked_halfmap1, input_map_2=masked_halfmap2)
        
        freq = frequency_array(fsc_curve, apix=apix)
        
        fsc_curves[i] = [freq, fsc_curve]

    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    ax1.grid(True)
    ax2 = ax1.twiny()
    if colors is None:
        colors = cm.rainbow(np.linspace(0,1,len(fsc_curves)))
        
    for i,fsc_profile in enumerate(fsc_curves.values()):
        freq = fsc_profile[0]
        fsc_curve = fsc_profile[1]
        
        ax1.plot(freq,fsc_curve, c=colors[i])
    ax1.plot(freq,np.ones(len(fsc_curve))*0.143,'k--')  

        
    if title is not None:
        ax1.set_title(title)
    ax2.set_xticks(ax1.get_xticks())
    ax2.set_xbound(ax1.get_xbound())
    ax2.set_xticklabels([round(1/x,1) for x in ax1.get_xticks()])
    if legend is None:
        legend = ['FSC curve','FSC=0.143','FSC=0.5']
        
    ax1.legend(legend,fontsize=font)
    ax1.set_xlabel(r'$1/d [\AA^{-1}]$',fontsize=font)
    ax1.set_ylabel('FSC',fontsize=font)
    ax2.set_xlabel('$d [\AA]$',fontsize=font)
    
    return fsc_curves
               
    

def plot_fsc_maps(input_map_1, input_map_2, apix, input_mask=None, calc_fsc=None,font=16, legend=None, title=None):

    import matplotlib.pyplot as plt
    from scipy.interpolate import interp1d
    from locscale.include.emmer.ndimage.map_utils import parse_input
    from locscale.include.emmer.ndimage.profile_tools import frequency_array
    
    emmap_1 = parse_input(input_map_1)
    emmap_2 = parse_input(input_map_2)
    
    if input_mask is not None:
        mask = parse_input(input_mask)
        masked_emmap_1 = mask*emmap_1
        masked_emmap_2 = mask*emmap_2
        fsc = calculate_fsc_maps(masked_emmap_1, masked_emmap_2)
    else:
        fsc = calculate_fsc_maps(emmap_1, emmap_2)
    freq = frequency_array(fsc, apix=apix)
    
    ## get relation between freq and fsc
    f = interp1d(fsc, freq)
    

    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    ax1.grid(True)
    ax2 = ax1.twiny()
    
    ax1.plot(freq,fsc, 'b')
    ax1.plot(freq,np.ones(len(fsc))*0.143,'k--')  
    ax1.plot(freq,np.ones(len(fsc))*0.5,'k-')  
    
    if title is not None:
        ax1.set_title(title)
    ax2.set_xticks(ax1.get_xticks())
    ax2.set_xbound(ax1.get_xbound())
    ax2.set_xticklabels([round(1/x,1) for x in ax1.get_xticks()])
    if legend is None:
        legend = ['FSC curve','FSC=0.143','FSC=0.5']
    
    ax1.legend(legend,fontsize=font)
    ax1.set_xlabel(r'$1/d [\AA^{-1}]$',fontsize=font)
    ax1.set_ylabel('FSC',fontsize=font)
    ax2.set_xlabel('$d [\AA]$',fontsize=font)
    
    
    if calc_fsc is not None:
        fsc_value = f(calc_fsc)
    
        return fig, fsc_value
    else:
        return fig

def plot_phase_correlation_maps(input_map_1, input_map_2, apix, input_mask=None, calc_fsc=None,font=16, legend=None, title=None):

    import matplotlib.pyplot as plt
    from locscale.include.emmer.ndimage.map_utils import parse_input
    from locscale.include.emmer.ndimage.profile_tools import frequency_array
    
    emmap_1 = parse_input(input_map_1)
    emmap_2 = parse_input(input_map_2)
    
    if input_mask is not None:
        mask = parse_input(input_mask)
        masked_emmap_1 = mask*emmap_1
        masked_emmap_2 = mask*emmap_2
        fsc = calculate_phase_correlation_maps(masked_emmap_1, masked_emmap_2)
    else:
        fsc = calculate_phase_correlation_maps(emmap_1, emmap_2)
    freq = frequency_array(fsc, apix=apix)
    


    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    ax1.grid(True)
    ax2 = ax1.twiny()
    
    ax1.plot(freq,fsc, 'b')

    
    if title is not None:
        ax1.set_title(title)
    ax2.set_xticks(ax1.get_xticks())
    ax2.set_xbound(ax1.get_xbound())
    ax2.set_xticklabels([round(1/x,1) for x in ax1.get_xticks()])
    if legend is None:
        legend = ['Phase correlation curve']
    
    ax1.legend(legend,fontsize=font)
    ax1.set_xlabel(r'$1/d [\AA^{-1}]$',fontsize=font)
    ax1.set_ylabel('$\Delta\phi$ ($1/d$)',fontsize=font)
    ax2.set_xlabel('$d [\AA]$',fontsize=font)

    return fig

def plot_amplitude_correlation_maps(input_map_1, input_map_2, apix, input_mask=None, calc_fsc=None,font=16, legend=None, title=None):

    import matplotlib.pyplot as plt
    from locscale.include.emmer.ndimage.map_utils import parse_input
    from locscale.include.emmer.ndimage.profile_tools import frequency_array
    
    emmap_1 = parse_input(input_map_1)
    emmap_2 = parse_input(input_map_2)
    
    if input_mask is not None:
        mask = parse_input(input_mask)
        masked_emmap_1 = mask*emmap_1
        masked_emmap_2 = mask*emmap_2
        fsc = calculate_amplitude_correlation_maps(masked_emmap_1, masked_emmap_2)
    else:
        fsc = calculate_amplitude_correlation_maps(emmap_1, emmap_2)
    freq = frequency_array(fsc, apix=apix)
    


    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    ax1.grid(True)
    ax2 = ax1.twiny()
    
    ax1.plot(freq,fsc, 'b')

    
    if title is not None:
        ax1.set_title(title)
    ax2.set_xticks(ax1.get_xticks())
    ax2.set_xbound(ax1.get_xbound())
    ax2.set_xticklabels([round(1/x,1) for x in ax1.get_xticks()])
    if legend is None:
        legend = ['Amplitude correlation curve']
    
    ax1.legend(legend,fontsize=font)
    ax1.set_xlabel(r'$1/d [\AA^{-1}]$',fontsize=font)
    ax1.set_ylabel('Amplitude',fontsize=font)
    ax2.set_xlabel('$d [\AA]$',fontsize=font)

    return fig
    
    
    


def get_fsc_filter(input_map_1, input_map_2):
    fsc_curve = calculate_fsc_maps(input_map_1, input_map_2)
    C_ref = 2*fsc_curve / (1+fsc_curve)
    
    return C_ref
    
def apply_fsc_filter(emmap, apix, fsc_curve=None, halfmap_1=None, halfmap_2=None):
    from locscale.include.emmer.ndimage.profile_tools import compute_radial_profile
    from locscale.include.emmer.ndimage.map_tools import compute_scale_factors, set_radial_profile
    if fsc_curve is not None:
        C_ref = 2*fsc_curve / (1+fsc_curve)
    elif halfmap_1 is not None and halfmap_2 is not None:
        C_ref = get_fsc_filter(halfmap_1, halfmap_2)
    
    
    rp_emmap, radii = compute_radial_profile(emmap, return_indices=True)
    fsc_filtered_rp = rp_emmap * C_ref
    
    sf = compute_scale_factors(rp_emmap, fsc_filtered_rp)
    filtered_emmap = set_radial_profile(emmap, sf, radii)
    
    return filtered_emmap, C_ref
    