#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 22 15:23:14 2021

@author: alok
"""

#from emmer.headers import *
import numpy as np


def frequency_array(amplitudes=None,apix=None,profile_size=None):
    '''
    Returns a numpy array with elements corresponding to the frequencies of a signal

    Parameters
    ----------
    amplitudes : numpy.ndarray (1,N)
        Amplitudes 
    apix : float
        pixel size, or more generally the size in real units for each index (time, or space)

    Returns
    -------
    freq : numpy.ndarray (1,N)
        Frequencies corresponding to the amplitudes, given the pixelsize
        

    '''
    if amplitudes is not None:
        n = len(amplitudes)
    elif profile_size is not None:
        n = profile_size
    else:
        print("Please enter the size of the array or send the array itself!")
        return 0
    
    if apix is None:
        print("Warning: voxelsize parameter not entered. \n Using apix = 1")
        apix = 1
        
    freq = np.linspace(1/(apix*n),1/(apix*2),n,endpoint=True)
    return freq

def plot_radial_profile(freq,list_of_profiles,legends=None, font=28,showlegend=True, showPoints=True, alpha=0.05, variation=None, logScale=True, ylims=None, xlims=None, crop_freq=None):
    import matplotlib.pyplot as plt
    from matplotlib.pyplot import cm
    from locscale.include.emmer.ndimage.profile_tools import crop_profile_between_frequency
    '''
    Given a list of amplitudes, plot them against a common frequency axis.

    Parameters
    ----------
    freq : np.ndarray
        Common frequency axis. Same size as the profiles in list of profiles
    list_of_profiles : list 
        List of amplitude profiles all having same size.
        list_of_profiles = [profile_1(type=np.ndarray), profile_2(type=np.ndarray), ...]
    colors : list, optional
        Custom color list. Max 6 entries. The default is ['r','g','b','k','y','m'].
    legends : lsit of string, optional
        Attach a legend corresponding to each profile in the list of profile. 
    font : int, optional
        fontsize for the plots. The default is 12.
    showlegend : bool, optional
        If you need to hide the legends, set this parameter to False

    Returns
    -------
    None.

    '''
    
    plt.rc('font',size=font)
        
    i = 0
    colors = cm.rainbow(np.linspace(0,1,len(list_of_profiles)))
    
   # if showPoints:
   #     colors = [x+".-" for x in colors]
   # else:
   #     colors = [x+"-" for x in colors]
    if legends is None:
        legends = ["profile_"+str(i) for i in range(len(list_of_profiles))]
    if len(list_of_profiles) <= 50:
        fig = plt.figure()
        ax1 = fig.add_subplot(111)
        ax1.grid(False)
        ax2 = ax1.twiny()
        #plt.xticks(fontsize=font)
        #plt.yticks(fontsize=font)
        if logScale:
            for profile in list_of_profiles:
                if crop_freq is not None:
                    freq, profile = crop_profile_between_frequency(freq, profile, crop_freq[0], crop_freq[1])
                ax1.plot(freq**2,np.log(profile),c=colors[i], linewidth=1)
                i += 1
            
            ax2.set_xticks(ax1.get_xticks())
            ax2.set_xbound(ax1.get_xbound())
            ax2.set_xticklabels([round(1/np.sqrt(x),1) for x in ax1.get_xticks()])
            if showlegend:
                ax1.legend(legends)
            ax1.set_xlabel(r'Spatial Frequency, $d^{-2} (\AA^{-2})$')
            ax1.set_ylabel(r'$ln  \langle \mid F \mid \rangle $ ')
            ax2.set_xlabel(r'$d (\AA)$')
        else:
            for profile in list_of_profiles:
                if crop_freq is not None:
                    freq, profile = crop_profile_between_frequency(freq, profile, crop_freq[0], crop_freq[1])
                ax1.plot(freq,profile,c=colors[i], linewidth=1)
                i += 1
            
            
            if showlegend:
                ax1.legend(legends)
        
                
            ax2.set_xticks(ax1.get_xticks())
            ax2.set_xbound(ax1.get_xbound())
            ax2.set_xticklabels([round(1/x,1) for x in ax1.get_xticks()])
            
    
            ax1.set_xlabel(r'Spatial Frequency, $1/d [\AA^{-1}]$')
            ax1.set_ylabel(r'Normalised $ \langle F \rangle $')
            ax2.set_xlabel('$d [\AA]$')
            
        
        
    elif variation is None:
        
        profile_list = np.array(list_of_profiles)
        average_profile = np.einsum("ij->j", profile_list) / len(profile_list)
        
        variation = []
        for col_index in range(profile_list.shape[1]):
            col_extract = profile_list[:,col_index]
            variation.append(col_extract.std())

        variation = np.array(variation)
        
        y_max = average_profile + variation
        y_min = average_profile - variation

        fig = plt.figure()
        
        ax1 = fig.add_subplot(111)
        ax1.grid(False)
        ax2 = ax1.twiny()
        
        if logScale:
            if crop_freq is not None:
                freq, average_profile = crop_profile_between_frequency(freq, average_profile, crop_freq[0], crop_freq[1])
                freq, y_max = crop_profile_between_frequency(freq, y_max, crop_freq[0], crop_freq[1])
                freq, y_min = crop_profile_between_frequency(freq, y_min, crop_freq[0], crop_freq[1])
            
            ax1.plot(freq**2, np.log(average_profile), 'k',alpha=1)
            ax1.fill_between(freq**2,np.log(y_max), np.log(y_min), color="grey", alpha=0.5)
            if showlegend:
                ax1.legend(["N={}".format(len(profile_list))])
        
            
            ax2.set_xticks(ax1.get_xticks())
            ax2.set_xbound(ax1.get_xbound())
            ax2.set_xticklabels([round(1/np.sqrt(x),1) for x in ax1.get_xticks()])
            
    
            ax1.set_xlabel(r'Spatial Frequency $1/d^2 [\AA^{-2}]$',fontsize=font)
            ax1.set_ylabel(r'$ln\mid F \mid $',fontsize=font)
            ax2.set_xlabel(r'$d [\AA]$',fontsize=font)
        else:
            ax1.plot(freq, average_profile, 'k',alpha=1)
            ax1.fill_between(freq,y_max, y_min,color="grey", alpha=0.5)
            
            if showlegend:
                ax1.legend(["N={}".format(len(profile_list))])
        
                
            ax2.set_xticks(ax1.get_xticks())
            ax2.set_xbound(ax1.get_xbound())
            ax2.set_xticklabels([round(1/x,1) for x in ax1.get_xticks()])
            
    
            ax1.set_xlabel(r'Spatial Frequency $1/d [\AA^{-1}]$',fontsize=font)
            ax1.set_ylabel(r'normalised $ \langle F \rangle $',fontsize=font)
            ax2.set_xlabel(r'$d [\AA]$',fontsize=font)
        
    else:
        if variation is None:
            raise UserWarning("Include a variation variable to plot radial profile with variance")
        
        if len(list_of_profiles) > 1:
            raise UserWarning("Multiple profiles given as average profile..")
            
        average_profile = list_of_profiles[0]
        
        y_max = average_profile + variation
        y_min = average_profile - variation

        fig = plt.figure()
        
        ax1 = fig.add_subplot(111)
        ax1.grid(True)
        ax2 = ax1.twiny()
        
        ax1.plot(freq**2, np.log(average_profile), 'k',alpha=1)
        ax1.fill_between(freq**2,np.log(y_max), np.log(y_min),color="grey", alpha=0.5)
        ax1.legend(["N={}".format(len(profile_list))])
        
        ax1.tick_params(axis='both',which='major')
        ax2.tick_params(axis='both',which='major')
        ax2.set_xticks(ax1.get_xticks())
        ax2.set_xbound(ax1.get_xbound())
        ax2.set_xticklabels([round(1/np.sqrt(x),1) for x in ax1.get_xticks()])
        

        ax1.set_xlabel(r'$1/d^2  [\AA^{-2}]$',fontsize=font)
        ax1.set_ylabel('$ln\mid F \mid $',fontsize=font)
        ax2.set_xlabel('$d [\AA]$',fontsize=font)

    
    if ylims is not None:
        plt.ylim(ylims)
        plt.yticks([-10,-5, 0])
    if xlims is not None:
        plt.xlim(xlims)
    
    
    plt.tight_layout()
    return fig
def plot_emmap_section(emmap, title="EMMAP Sections"):
    import matplotlib.pyplot as plt
    
    fig = plt.figure()
    plt.title(title)
    plt.axis("off")
    zn, yn, xn = emmap.shape
    
    ax1 = fig.add_subplot(131)
    plt.title("XY plane")
    plt.imshow(emmap[zn//2,:,:])
    
    ax2 = fig.add_subplot(132)
    plt.title("YZ plane")
    plt.imshow(emmap[:,:,xn//2])
    ax2.yaxis.set_major_locator(plt.NullLocator())

    
    ax3 = fig.add_subplot(133)
    plt.title("XZ plane")
    plt.imshow(emmap[:,yn//2,:])
    ax3.yaxis.set_major_locator(plt.NullLocator())

    
    return fig
    
    
def add_deviations_to_reference_profile(freq, reference_profile, scaled_theoretical_profile, wilson_cutoff, nyquist_cutoff, deviation_freq_start, deviation_freq_end, magnify=1):
    '''
    Function to add deviations from a reference profile which is assumed to be exponential at high frequencies

    Parameters
    ----------
    freq : TYPE
        DESCRIPTION.
    reference_profile : TYPE
        DESCRIPTION.
    scaled_theoretical_profile : TYPE
        DESCRIPTION.
    deviation_freq_start : TYPE
        DESCRIPTION.
    deviation_freq_end : TYPE
        DESCRIPTION.

    Returns
    -------
    deviated_profile_tuple : tuple
        (freq, deviated_profile)

    '''
    
    deviations, exponential_fit = calculate_required_deviation(freq, scaled_theoretical_profile, wilson_cutoff, nyquist_cutoff, deviation_freq_start, deviation_freq_end)
    if magnify > 1:
        #f = magnification_function(magnify)
        #deviated_reference_profile = reference_profile * f(deviations)
        deviated_reference_profile = reference_profile + magnify * deviations
    else:
        #deviated_reference_profile = np.log(reference_profile) * deviations
        deviated_reference_profile = reference_profile + deviations
    
    return deviated_reference_profile, exponential_fit
    
def magnification_function(magnify, cutoff=1, x_max = 10):
    '''
    Returns a magnified deviations curve based on product 

    Parameters
    ----------
    deviations : TYPE
        DESCRIPTION.
    magnify : TYPE
        DESCRIPTION.
    cutoff : TYPE, optional
        DESCRIPTION. The default is 1.

    Returns
    -------
    None.

    '''
    from scipy.interpolate import interp1d
    xdata = np.array(list(np.linspace(0, cutoff, 100))+list(np.linspace(cutoff, x_max,100)))
    ydata = []
    for x in xdata:
        if x < cutoff:
            ydata.append(x/magnify)
        elif x > cutoff:
            ydata.append(x*magnify)
        else:
            ydata.append(cutoff)
    
    ydata = np.array(ydata)
    
    f = interp1d(x=xdata, y=ydata)
    
    return f
    
def merge_two_profiles(profile_1,profile_2,freq, smooth=1, d_cutoff=None, f_cutoff=None):
    '''
    Function to merge two profiles at a cutoff threshold based on differential weighting of two profiles

    Parameters
    ----------
    profile_1 : numpy.ndarray
        
    profile_2 : numpy.ndarray
        same size of profile_1
    freq : numpy.ndarray
        Frequencies corresponding to both profile_1 and profile_2
    d_cutoff : float
        Cutoff frequency defined in terms of distance (unit = A)
    f_cutoff : float
        Cutoff frequency given in terms of spatial frequency (unit = 1/A)
    smooth : float, optional
        smoothening parameter to control the transition region of two profiles

    Returns
    -------
    merged_profile : tuple of two numpy.ndarray
    
    '''

    if not (len(freq) == len(profile_1) and len(freq) == len(profile_2)):
        print("Size of two profiles not equivalent. Please check the dimensions and give another input")
        return None
    
    k = smooth
    d = 1 / freq
    
    if d_cutoff is not None:
        d_cutoff = d_cutoff
    
    elif f_cutoff is not None:
        d_cutoff = 1 / f_cutoff
    
    else:
        print("Please enter a cutoff frequency either in terms of spatial frequency (1/A) or distance (A)")
        return None
    
    weight_1 = 1 / (1 + np.exp(k * (d_cutoff - d)))
    weight_2 = 1 - weight_1
    
    merged_profile = weight_1 * profile_1 + weight_2 * profile_2
    
    return merged_profile

def compute_radial_profile(vol, center=[0,0,0], return_indices=False):
    '''
    Computes the radial profile of a given volume

    Parameters
    ----------
    vol : numpy.ndarray
        Input array
    center : list, optional
        DESCRIPTION. The default is [0,0,0].
    return_indices : bool, optional
        

    Returns
    -------
    radial_profile : numpy.ndarray (1D)
        Radial profile
        

    '''
    dim = vol.shape
    m = np.mod(vol.shape,2)
    # make compliant with both fftn and rfftn
    if center is None:
        ps = np.abs(np.fft.fftshift((np.fft.fftn(vol))))
        z, y, x = np.indices(ps.shape)
        center = tuple((a - 1) / 2.0 for a in ps.shape[::-1])
        radii = np.sqrt((x - center[0])**2 + (y - center[1])**2 + (z - center[2])**2)
        radii = radii.astype(np.int)
    else:
        ps = np.abs( np.fft.rfftn(vol) )
        if not return_indices:
            x, y, z = np.indices(ps.shape)
            radii = np.sqrt(x**2 + y**2 + z**2)
            radii = radii.astype(np.int)
        else:
            [x, y, z] = np.mgrid[-dim[0]//2+m[0]:(dim[0]-1)//2+1, -dim[1]//2+m[1]:(dim[1]-1)//2+1, 0:dim[2]//2+1]
            x = np.fft.ifftshift(x)
            y = np.fft.ifftshift(y)
            radii = np.sqrt(x**2 + y**2 + z**2)
            radii = radii.astype(np.int)
    radial_profile = np.bincount(radii.ravel(), ps.ravel()) / np.bincount(radii.ravel())
    # exclude corner frequencies
    radial_profile = radial_profile[0:int(ps.shape[0]/2)]
    if not return_indices:
        return radial_profile
    else:
        return radial_profile, radii

def offset_radial_profile(vol, offset, radii):
    ps = np.fft.rfftn(vol)
    for j,r in enumerate(np.unique(radii)[0:vol.shape[0]//2]):
            idx = radii == r
            ps[idx] += offset

    return np.fft.irfftn(ps, s=vol.shape)
    
def compute_radial_profile_from_mrcs(mrc_paths,keys=None,logScale=False, ylims=None, xlims=None, crop_freq=None):
    '''
    Given a list of mrc paths, this function will extract volume data and plots a radial profile for each map. The legend by default will be the filename of the MRC map. Max six maps will be used as inputs (limit from plot_radial_profiles function)

    Parameters
    ----------
    mrc_paths : list of strings
        ["path/to/map1.mrc", "path/to/map2.mrc",..]
    keys : list of strings, optional
        ["map A", "map B",..]. The default is None.

    Returns
    -------
    radial_profiles : dict 
        Dictionary of radial profiles for each map.
    emmaps : dict
        Dictionary of emmap volumes for each map.

    '''
    import mrcfile
    
    if keys is None:
        keys = [path.split('/')[-1] for path in mrc_paths]

    mrcs = []
    for mrc in mrc_paths:
        mrcs.append(mrcfile.open(mrc))
    
    k = 0
    emmaps = {}
    radial_profiles = {}
    freq={}
    
    for mrc in mrcs:
        emmaps[keys[k]] = mrc.data
        k += 1
    for key in keys:
        radial_profiles[key] = compute_radial_profile(emmaps[key])
        radial_profiles[key] = radial_profiles[key]/radial_profiles[key].max()
        
    k = 0
    for key in keys:
        shapes = radial_profiles[key].shape[0]
        apix = mrcs[k].voxel_size.x
        freq[key] = np.linspace(1./(float(apix)*shapes), 1./(float(apix)*2), shapes,endpoint=True)
        k += 1 
    
    
    fig=plot_radial_profile(freq[keys[0]], list(radial_profiles.values()),legends=keys, logScale=logScale, showPoints=False, ylims=ylims, crop_freq=crop_freq,  xlims=xlims)
    
    for key in keys:
        radial_profiles[key] = tuple([freq,radial_profiles[key]])
        
    return fig
            

def measure_debye(freqs,amplitudes):
    '''
    Function to measure the "Debye Effect" from a radial profile

    Parameters
    ----------
    freqs : numpy.ndarray
        Frequency array
    amplitudes : numpy.ndarray
        Amplitudes

    Returns
    -------
    Debye effect : dict 
        Dictionary of debye effects at different identified peaks
    freq_step : float
        
    filtered amplitudes : numpy.ndarray
        DESCRIPTION.
    exponential fit : numpy.ndarray
        DESCRIPTION.

    '''
    from scipy.optimize import curve_fit
    from scipy import signal
    from emmer.ndimage.filter import butter_lowpass_filter
    
    # First filter the amplitude profile using Butterworth filter at nyq/2 and order =1
    n = len(amplitudes)
    freq_step = (freqs[-1]-freqs[0])/n
    
    nyq_freq = 1/(freq_step * 2)
    
    filtered_amplitudes = butter_lowpass_filter(amplitudes, nyq_freq/2, nyq_freq,order=1)
    
    debye_peaks = signal.find_peaks(filtered_amplitudes)[0]
    
    #debye_freqs = debye_peaks * freq_step
    
    
    ## Find the amplitude of the spectrum at debye frequency
    
    # First get a best match exponential fit
    
    def exponential(x,a,b):
        return a*np.exp(x*b)
    
    optimised_parameters,covariance = curve_fit(exponential,freqs,filtered_amplitudes)
    fit_a, fit_b = optimised_parameters
    #print(optimised_parameters)
    exponential_fit_amplitudes = exponential(freqs, fit_a, fit_b)
    debye_effect = {}
    for peak_index in debye_peaks:
        if peak_index*freq_step < 0.4:
            exponential_amplitude_at_debye_freq = exponential_fit_amplitudes[peak_index]
            filtered_amplitude_at_debye_freq = filtered_amplitudes[peak_index]
            debye_effect[peak_index] = filtered_amplitude_at_debye_freq - exponential_amplitude_at_debye_freq
    
    return debye_effect,freq_step,filtered_amplitudes,exponential_fit_amplitudes

def estimate_bfactor_standard(freq, amplitude, wilson_cutoff, fsc_cutoff, return_amplitude=False, return_fit_quality=False, standard_notation=False):
    '''
    From a given radial profile, estimate the b_factor from the high frequency cutoff

    Parameters
    ----------
    freq : numpy.ndarray
        Frequency array
    amplitude : numpy.ndarray
        Amplitudes
    wilson_cutoff : float
        Frequency from which wilson statistics are valid. Units: Angstorm
    fsc_cutoff : float
        FSC resolution calculated at 0.143 (for halfmaps). Units: Angstorm
        

    Returns
    -------
    b_factor : float
        The estimated b factor
    
    amp : float
        The estimated amplitude of the exponential fit

    '''
    from scipy.optimize import curve_fit
    from sklearn.metrics import r2_score
    
    def linear_fit(xdata,slope,const):
        ydata = const + slope*xdata
        return ydata
    
    wilson_freq = 1 / wilson_cutoff
    fsc_freq = 1 / fsc_cutoff
    
    if freq[0] >= wilson_freq:
        start_index = 0
    else:
        start_index = np.where(freq>=wilson_freq)[0][0]
    
    if freq[-1] <= fsc_freq:
        end_index = len(freq)
    else:
        end_index = np.where(freq>=fsc_freq)[0][0]
    
    xdata = freq[start_index:end_index]**2
    ydata = np.log(amplitude[start_index:end_index])
    
    
    param, _ = curve_fit(linear_fit,xdata,ydata)
    
    if standard_notation:
        b_factor = -1 * param[0] * 4   ## Inverse of slope
    else:
        b_factor = param[0] * 4
    
    exp_fit_amplitude = np.exp(param[1])
    
    #print("B factor: "+str(round(param[0]*4,2)))
    y_pred = linear_fit(xdata, slope=param[0], const=param[1])
    r2 = r2_score(y_true=ydata, y_pred=y_pred)
    if return_amplitude:
        if return_fit_quality:
            return b_factor,exp_fit_amplitude, r2
        else:
            return b_factor,exp_fit_amplitude
    else:
        if return_fit_quality:
            return b_factor, r2
        else:
            return b_factor

def calculate_required_deviation(freq, scaled_theoretical_profile, wilson_cutoff, nyquist_cutoff, deviation_freq_start, deviation_freq_end=None):
    '''
    Function to calculate the deviations per frequency from a scaled theoretical profile

    Parameters
    ----------
    scaled_theoretical_profile : TYPE
        DESCRIPTION.
    wilson_cutoff : TYPE
        DESCRIPTION.
    fsc_cutoff : TYPE
        DESCRIPTION.

    Returns
    -------
    deviations : numpy.ndarray
    deviations = scaled_theoretical - exponential_fit

    '''
    bfactor, amp = estimate_bfactor_standard(freq, scaled_theoretical_profile, wilson_cutoff=wilson_cutoff, fsc_cutoff=nyquist_cutoff, 
                                             return_amplitude=True)
    
    exponential_fit = amp * np.exp(bfactor * 0.25 * freq**2)
    
    #deviations = scaled_theoretical_profile / exponential_fit
    deviations = scaled_theoretical_profile - exponential_fit
    
    deviation_freq_start_freq = 1/deviation_freq_start
    
    start_index = np.where(freq>=deviation_freq_start_freq)[0][0]
    #deviations[:start_index] = 1
    deviations[:start_index] = 0
    if deviation_freq_end is not None:
        deviation_freq_end_freq = 1/deviation_freq_end
        end_index = np.where(freq>=deviation_freq_end_freq)[0][0]
        #deviations[end_index:] = 1
        deviations[end_index:] = 0
    
    return deviations, exponential_fit
    
    
    

def scale_profiles(reference_profile_tuple, scale_profile_tuple, wilson_cutoff, fsc_cutoff, return_bfactor_properties=False):
    '''
    Function to scale an input theoretical profile to a reference profile

    Parameters
    ----------
    reference_profile_tuple : tuple
        (freq_reference, amplitude_reference)
    scale_profile_tuple : tuple
        (freq_theoretical, amplitude_theoretical)
    just_use_exponential : bool, optional
        Returns just an exponential fit and not a scaled profile
    using_reference_profile : TYPE, optional
        DESCRIPTION. The default is False.
    start_freq : TYPE, optional
        DESCRIPTION. The default is 0.3.

    Returns
    -------
    TYPE
        DESCRIPTION.

    '''
    #power_zero_freq = reference_profile_tuple[1].max() 
    
    freq = reference_profile_tuple[0]
    reference_amplitude = reference_profile_tuple[1]
    
    freq_scale = scale_profile_tuple[0]
    scale_amplitude = scale_profile_tuple[1]
 
    bfactor_reference, fit_amp_reference, quality_of_fit = estimate_bfactor_standard(freq, reference_amplitude, wilson_cutoff=wilson_cutoff, fsc_cutoff=fsc_cutoff, return_amplitude=True,return_fit_quality=True)
    bfactor_scale, fit_amp_scale = estimate_bfactor_standard(freq_scale, scale_amplitude, wilson_cutoff=wilson_cutoff, fsc_cutoff=fsc_cutoff, return_amplitude=True)
    
    bfactor_diff = bfactor_reference-bfactor_scale
    
    amp_scaling_factor = fit_amp_reference / fit_amp_scale
        
    amplitude_scaled = amp_scaling_factor * scale_amplitude * np.exp(bfactor_diff * freq**2 / 4)
    
    
    if return_bfactor_properties:
        return (freq,amplitude_scaled), (bfactor_reference, fit_amp_reference, quality_of_fit)
    else:
        return (freq, amplitude_scaled)
    
def get_debye_statistics(radial_profiles_using_pdb): 
    import sys
    
    helix_freqs = []
    sheet_freqs = []
    missed = 0
    min_debye_effect_mag = 0.005
    max_debye_effect_mag = 0.01
    
    min_debye_freq = 0.1
    max_debye_freq = 0.3
    helix_profile_info = {}
    sheet_profile_info = {}
       
    for pdb in radial_profiles_using_pdb.keys():
        
        try:
            #apix = float(radial_profiles_using_pdb[pdb]['apix'])
            #factor = 0.25/apix
            
            #freq = np.array(radial_profiles_using_pdb[pdb]['freq'])
            #print(freq[0],freq[-1],len(freq))
            
            helix_amplitude=np.array(radial_profiles_using_pdb[pdb]['helix']) / max(radial_profiles_using_pdb[pdb]['helix'])
            sheet_amplitude=np.array(radial_profiles_using_pdb[pdb]['sheet']) / max(radial_profiles_using_pdb[pdb]['sheet'])
            
            ## This for new data format in radial_profiles_using_pdb starting from March 2021 
            apix_array = radial_profiles_using_pdb[pdb]['apix'][0]
            #print(apix_array.x)
            if apix_array.x == apix_array.y and apix_array.x == apix_array.z:
                #print("Same")
                apix = apix_array.x
            else:
                print("different for "+pdb)
            n = len(helix_amplitude)
            freq = np.linspace(1/(apix*n),1/(apix*2),n,endpoint=True)
            
            
            debye_effect_helix,freq_step_helix,_,_ = measure_debye(freq,helix_amplitude)
            debye_effect_sheet,freq_step_sheet,_,_ = measure_debye(freq,sheet_amplitude)
            
            for possible_debye_freq in debye_effect_helix.keys():
                if min_debye_freq <= possible_debye_freq*freq_step_helix <= max_debye_freq and min_debye_effect_mag <= debye_effect_helix[possible_debye_freq] <= max_debye_effect_mag and helix_amplitude[possible_debye_freq] < 0.02:
                    helix_freqs.append(possible_debye_freq*freq_step_helix)
                    helix_profile_info[pdb] = tuple([[freq,radial_profiles_using_pdb[pdb]['helix']],debye_effect_helix])
            for possible_debye_freq in debye_effect_sheet.keys():
                if min_debye_freq <= possible_debye_freq*freq_step_sheet <= max_debye_freq and min_debye_effect_mag <= debye_effect_sheet[possible_debye_freq] <= max_debye_effect_mag and sheet_amplitude[possible_debye_freq] < 0.02:
                    sheet_freqs.append(possible_debye_freq*freq_step_sheet)        
                    sheet_profile_info[pdb] = tuple([[freq,radial_profiles_using_pdb[pdb]['sheet']],debye_effect_sheet])
        except:
            missed += 1
            print(pdb,sys.exc_info())
            pass
            
    return [helix_freqs,sheet_freqs,helix_profile_info,sheet_profile_info,missed]
        
def resample_1d(x_old,y_old,num,xlims=None):
    '''
    Sample an given x-y data in a new grid 

    Parameters
    ----------
    x_old : numpy.ndarray
        data in x axis (same dim as y_old)
    y_old : numpy.ndarray
        data in y axis (same dim as x_old)
    num : int
        new number of data points

    Returns
    -------
    x_new : numpy.ndarray
        resampled x axis
    y_new : numpy.ndarray
        resampled y axis

    '''
    from scipy.interpolate import interp1d

    f = interp1d(x_old, y_old,kind='slinear',fill_value='extrapolate')
    if xlims is None:
        x_new = np.linspace(x_old[0], x_old[-1], num=num)
    else:
        xmin = xlims[0]
        xmax = xlims[1]
        x_new = np.linspace(xmin, xmax,num=num)
        
    y_new = f(x_new)
    #y_new[y_new>1]=1
    return x_new, y_new
    
def average_profiles(profiles_dictionary,num=1000):
    from locscale.include.emmer.ndimage.filter import get_nyquist_limit, butter_lowpass_filter, fit_series
    import sys
    
    resampled_profiles = {}
    min_freq,max_freq = find_xmin_xmax([data[0] for data in profiles_dictionary.values()])
    
    for pdb in profiles_dictionary.keys():
        try:
            
            freq_h = profiles_dictionary[pdb][0][0]
            nyq = get_nyquist_limit(freq_h)
            amplitudes_nofit_nonorm = butter_lowpass_filter(profiles_dictionary[pdb][0][1],nyq/2,nyq,1)
                  
            freq_h,amp_fit_nonorm = fit_series([freq_h,amplitudes_nofit_nonorm],min_freq,max_freq,num)     
            amp_fit_norm = amp_fit_nonorm / amp_fit_nonorm.max()
            
            #resampled_freq,resampled_amp = resample_1d(freq_h,amp_fit_norm, num, xlims=[min_freq,max_freq])
            resampled_profiles[pdb] = [freq_h,amp_fit_norm]
            #plt.plot(freq_h,amp_fit_norm,'k'), 
            
        except:
            
            e = sys.exc_info()[0]
            f = sys.exc_info()[1]
            o = sys.exc_info()[2]
            print(pdb,e,f,o)
    
    average_profile = np.zeros(num)
    for pdb in resampled_profiles.keys():
        average_profile += resampled_profiles[pdb][1]
    common_freq = np.linspace(min_freq, max_freq,num)
    average_profile /= len(profiles_dictionary.keys())
    
    return common_freq, average_profile, resampled_profiles    
    
   
def find_xmin_xmax(profiles):
    '''
    profiles: python.list containing profile
    
    profile = [x_data,y_data]
    profiles = [profile1,profile2,profile3...]
    
    '''
    
    for i,profile in enumerate(profiles):
        
        if i == 0:
            xmin = profile[0][0]
            xmax = profile[0][-1]
            
        else:
            if profile[0][0] < xmin:
                xmin = profile[0][0]
            if profile[0][-1] > xmax:
                xmax = profile[0][-1]
    return xmin,xmax    
    
def number_of_segments(fsc_resolution):
    if fsc_resolution < 3:
        return 4
    elif fsc_resolution >= 3 and fsc_resolution < 5:
        return 3
    elif fsc_resolution >= 5 and fsc_resolution < 6:
        return 2
    else:
        print("Warning: resolution too low to estimate cutoffs. Returning 1")
        return 1

def plot_pwlf_fit(emmap_path, mask_path, fsc_resolution):
    '''
    Function to plot PWLF fit for a given input map

    Parameters
    ----------
    emmap_path : TYPE
        DESCRIPTION.
    mask_path : TYPE
        DESCRIPTION.
    fsc_resolution : TYPE
        DESCRIPTION.

    Returns
    -------
    None.

    '''    
    import mrcfile
    import os
    from locscale.include.emmer.ndimage.profile_tools import compute_radial_profile, frequency_array, plot_radial_profile
    from locscale.include.emmer.pdb.pdb_tools import find_wilson_cutoff
    
    emmap_name = os.path.basename(emmap_path)
    
    emmap = mrcfile.open(emmap_path).data
    apix = mrcfile.open(emmap_path).voxel_size.tolist()[0]
    rp_emmap = compute_radial_profile(emmap)
    
    freq = frequency_array(rp_emmap, apix=apix)
    wilson_cutoff = find_wilson_cutoff(mask_path=mask_path)
    
    bfactor, amplitude_zero_freq, (piecewise_linfit, z, slopes) = estimate_bfactor_through_pwlf(freq, rp_emmap, wilson_cutoff, fsc_cutoff=fsc_resolution, return_all=True)
    
    rp_fit = np.exp(piecewise_linfit.predict(freq**2))  ## Fit was trained using log scale data
    
    fig = plot_radial_profile(freq, [rp_emmap, rp_fit], legends=[emmap_name, "PWLF prediction"])
    print("Breakpoints at: {}".format((1/np.sqrt(z)).round(2)))
    
    return fig
    

def crop_profile_between_frequency(freq, amplitude, start_cutoff, end_cutoff):
    start_freq = 1 / start_cutoff
    end_freq = 1/end_cutoff
        
    if freq[0] >= start_freq:
        start_index = 0
    else:
        start_index = np.where(freq>=start_freq)[0][0]
        
    if freq[-1] <= end_freq:
        end_index = len(freq)
    else:
        end_index = np.where(freq>=end_freq)[0][0]
    
    crop_freq = freq[start_index:end_index]
    crop_amplitude = amplitude[start_index:end_index]
    
    return crop_freq, crop_amplitude
    
def estimate_bfactor_through_pwlf(freq,amplitudes,wilson_cutoff,fsc_cutoff, return_all=True, num_segments=None):
    '''
    Function to automatically find out linear region in a given radial profile 


    @Manual{pwlf,
            author = {Jekel, Charles F. and Venter,     Gerhard},
            title = {{pwlf:} A Python Library for Fitting 1D Continuous Piecewise Linear Functions},
            year = {2019},
            url = {https://github.com/cjekel/piecewise_linear_fit_py}
}

    Parameters
    ----------
    freq : numpy.ndarray
        
    amplitudes : numpy.ndarray
        

    Returns
    -------
    start_freq_in_angstorm, estimated_bfactor

    '''
    import pwlf
    from locscale.pseudomodel.pseudomodel_headers import number_of_segments
    
    if num_segments is None:
            num_segments = number_of_segments(fsc_cutoff)
            
    if num_segments < 2:
        print("Number of segments = 1 using standard method of evaluating bfactor")
        bfactor, amplitude_zero_freq = estimate_bfactor_standard(freq, amplitudes, wilson_cutoff, fsc_cutoff, return_amplitude=True)
        piecewise_linfit = amplitude_zero_freq * np.exp(0.25 * bfactor * freq**2)
        z = [(1/wilson_cutoff)**2, (1/fsc_cutoff)**2]
        slopes = [bfactor / 4]
    
    else:
        
        start_freq = 1 / wilson_cutoff
        end_freq = 1/fsc_cutoff
        
        if freq[0] >= start_freq:
            start_index = 0
        else:
            start_index = np.where(freq>=start_freq)[0][0]
        
        if freq[-1] <= end_freq:
            end_index = len(freq)
        else:
            end_index = np.where(freq>=end_freq)[0][0]
        
        x_data = freq[start_index:end_index]**2
        y_data = np.log(amplitudes[start_index:end_index])
        
        piecewise_linfit = pwlf.PiecewiseLinFit(x_data, y_data)
        
        z = piecewise_linfit.fit(n_segments=num_segments)
        
        slopes = piecewise_linfit.calc_slopes()
        
        bfactor = slopes[-1] * 4
        
        amplitude_zero_freq = piecewise_linfit.predict(0)
    
    if return_all:
        return bfactor, amplitude_zero_freq, (piecewise_linfit, z, slopes)
    else:
        return bfactor

def get_theoretical_profile(length,apix, profile_type='helix'):
    import pickle
    from locscale.include.emmer.ndimage.profile_tools import resample_1d
    from locscale.pseudomodel.pseudomodel_headers import check_dependencies
    import os
    
    path_to_locscale = check_dependencies()['locscale']
    location_of_theoretical_profiles = os.path.join(path_to_locscale, "locscale","utils","theoretical_profiles.pickle")
    
    with open(location_of_theoretical_profiles,'rb') as f:
        profiles = pickle.load(f)
    
    frequency_limits = (float(1/(apix*length)),float(1/(apix*2)))
    helix_profile = profiles[profile_type]
    resampled_helix_profile = resample_1d(helix_profile['freq'], helix_profile['profile'],num=length,xlims=frequency_limits)
    return resampled_helix_profile
    
    
    

def generate_no_debye_profile(freq, amplitudes, wilson_cutoff=10, smooth=1):
    from locscale.include.emmer.ndimage.profile_tools import merge_two_profiles, estimate_bfactor_standard
    
    bfactor, amp = estimate_bfactor_standard(freq, amplitudes, wilson_cutoff=10, fsc_cutoff=1/freq[-1], return_amplitude=True, standard_notation=True)
    
    exponential_fit = amp * np.exp(-0.25 * bfactor * freq**2)

    y_data = np.log(amplitudes)
    x_data = freq**2
    y_fit = np.log(exponential_fit)
    
    #merged_profile = np.concatenate((ydata[:start_index],y_data_wilson))
    merged_profile = merge_two_profiles(y_data, y_fit, freq, d_cutoff=wilson_cutoff, smooth=smooth)
    
    new_amplitudes = np.exp(merged_profile)
    
    return new_amplitudes

