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
        deviated_reference_profile = reference_profile + magnify * deviations
    else:
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
        radii = radii.astype(int)
    else:
        ps = np.abs( np.fft.rfftn(vol) )
        if not return_indices:
            x, y, z = np.indices(ps.shape)
            radii = np.sqrt(x**2 + y**2 + z**2)
            radii = radii.astype(int)
        else:
            [x, y, z] = np.mgrid[-dim[0]//2+m[0]:(dim[0]-1)//2+1, -dim[1]//2+m[1]:(dim[1]-1)//2+1, 0:dim[2]//2+1]
            x = np.fft.ifftshift(x)
            y = np.fft.ifftshift(y)
            radii = np.sqrt(x**2 + y**2 + z**2)
            radii = radii.astype(int)
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
    
    ## TBC
def amplitude_from_resolution(freq, amplitudes, probe_resolution, logScale = True):
    from scipy.interpolate import interp1d
    
    if probe_resolution < 0:
        raise UserWarning("Enter probe resolution > 0A")
    if logScale:
        xdata = freq**2
        ydata = np.log(amplitudes)
        f = interp1d(xdata, ydata)
        probe_freq = 1/probe_resolution
        x_probe = probe_freq**2
        y_probe = f(x_probe)
        
        amplitude_probe = np.exp(y_probe)
        
        return amplitude_probe
        
    else:
        xdata = freq
        ydata = amplitudes
        f = interp1d(xdata, ydata)
        probe_freq = 1/probe_resolution
        x_probe = probe_freq
        y_probe = f(x_probe)
        amplitude_probe = y_probe
        
        return amplitude_probe
        
def resolution_from_amplitude(freq, amplitudes, probe_amplitude, logScale = True, suppress_warning=False):
    from scipy.interpolate import interp1d
    
    if probe_amplitude < 0:
        raise UserWarning("Enter probe amplitude > 0A")
        
    if logScale:
        xdata = np.log(amplitudes)
        ydata = freq**2
        g = interp1d(xdata, ydata)
        x_probe = np.log(probe_amplitude)
        if x_probe < xdata.min():
            return 1/freq[-1]
        if x_probe > xdata.max():
            return 1/freq[0]
        y_probe = g(x_probe)
        
        freq_probe = np.sqrt(y_probe)
        if freq_probe < freq[-1]:
            probe_resolution = 1/freq_probe
            return probe_resolution
        else:
            if not suppress_warning:
                print("Warning: estimated resolution outside interpolation range. Returning Nyquist {:.2f}".format(1/freq_probe))
            return 1/freq[-1]
    else:
        xdata = amplitudes
        ydata = freq
        g = interp1d(xdata, ydata)
        
        x_probe = probe_amplitude
        y_probe = g(x_probe)
        freq_probe = y_probe
        if freq_probe < freq[-1]:
            probe_resolution = 1/freq_probe
            return probe_resolution  
        else:
            if not suppress_warning:
                print("Warning: estimated resolution outside interpolation range. Returning Nyquist")
            return 1/freq[-1]
        
    
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
    import warnings
    from locscale.include.emmer.ndimage.profile_tools import number_of_segments
    
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
        ## Ignore RunTimeWarning: invalid value encountered in log
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", category=RuntimeWarning)
            piecewise_linfit = pwlf.PiecewiseLinFit(x_data, y_data)
            z = piecewise_linfit.fit(n_segments=num_segments, disp=False)
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
    from locscale.utils.file_tools import get_locscale_path
    import os
    
    path_to_locscale = get_locscale_path()
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

