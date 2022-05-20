import numpy as np

import mrcfile

def calculate_fourier_frequencies(im, apix):
    """Return the image frequency for every voxel in Fourierspace
    for n-dimensional images
    """
    per_axis_freq = [np.fft.fftfreq(N) for N in im.shape[:-1]]
    per_axis_freq.append(np.fft.rfftfreq(im.shape[-1]))
    dims = np.meshgrid(*per_axis_freq, indexing='ij', sparse=True)
    fourier_frequencies = np.sqrt(np.sum([dim**2 for dim in dims]))
    fourier_frequencies_angstrom = fourier_frequencies / apix
    return fourier_frequencies_angstrom

def critical_exposure(fourier_frequency):
    """Returns exposure (e-/A**2) for which Fourier information of a given spatial frequency
    shrinks to 1/e. Based on (Grant,Grigorieff 2015). 
    """
    return 0.245*fourier_frequency**(-1.665)+2.81 #Grant,Grigorieff 2015

def dose_weight_per_fourier_frequency(fourier_frequencies, accumulated_dose):
    """Returns weight for each fourier pixel based on accumulated dose.
    """
    critical_exposures = critical_exposure(fourier_frequencies)
    dose_weight = np.exp(-accumulated_dose/(2*critical_exposures))
    return dose_weight

def doseweight_image(im, dose, apix):
    im_freq         = calculate_fourier_frequencies(im, apix)
    im_exp_filter   = dose_weight_per_fourier_frequency(im_freq, dose)
    im_fft          = np.fft.rfft2(im)
    im_fft_filtered = im_fft * im_exp_filter
    im_filtered     = np.fft.irfft2(im_fft_filtered)
    return im_filtered

def tanh_filter(im_freq, cutoff):
    """Returns filter coefficients for a hyperbolic tangent filter. 
    """
    cutoff_freq = 1/cutoff
    filter_fall_off = 0.1;
    filter_coefficients = 1.0 - (1.0 - 0.5*(np.tanh((np.pi*(im_freq+cutoff_freq)/(2*filter_fall_off*cutoff_freq))) - np.tanh((np.pi*(im_freq-cutoff_freq)/(2*filter_fall_off*cutoff_freq)))));
    return filter_coefficients;

def low_pass_filter(im, cutoff, apix):
    """
    Returns a low-pass filter image from a tanh filter.
    """
    im_freq     = calculate_fourier_frequencies(im, apix=apix)
    im_filter   = tanh_filter(im_freq, cutoff);
    im_fft      = np.fft.rfftn(im)
    im_fft_filtered = im_fft * im_filter
    im_filtered = np.fft.irfftn(im_fft_filtered)
    return im_filtered

def high_pass_filter(im, cutoff, apix):
    """
    Returns a high-pass filter image from a tanh filter.
    """
    im_freq     = calculate_fourier_frequencies(im, apix=apix)
    im_filter   = 1-tanh_filter(im_freq, cutoff);
    im_fft      = np.fft.rfftn(im)
    im_fft_filtered = im_fft * im_filter
    im_filtered = np.fft.irfftn(im_fft_filtered)
    return im_filtered

def band_pass_filter(im, cutoff_lims, apix):
    '''
    

    Parameters
    ----------
    im : TYPE
        DESCRIPTION.
    cutoff_lims : list
        DESCRIPTION.
    apix : TYPE
        DESCRIPTION.

    Returns
    -------
    None.

    '''
    low_pass_filter_cutoff = min(cutoff_lims)
    high_pass_filter_cutoff = max(cutoff_lims)
    
    im_low_pass = low_pass_filter(im, low_pass_filter_cutoff, apix)
    im_high_pass = high_pass_filter(im_low_pass, high_pass_filter_cutoff, apix)
    
    return im_high_pass

def apply_filter_to_map(emmap_path,dmin,output_filename=None):
    '''
    This function applies a low pass tanh filter to a MRC file given the map and filter cutoff as input

    Parameters
    ----------
    emmap_path : string
        path/to/emmap.mrc
    dmin : float
        low pass filter cutoff to apply tanh filter

    Returns
    -------
    filtered_map : string
        path/to/filtered.mrc. By default, it would have same name as emmap path, with "_filtered_dmin_A.mrc"
        
        

    '''    
    from locscale.include.emmer.ndimage.map_utils import save_as_mrc
    from locscale.include.emmer.ndimage.filter import low_pass_filter
        
    emmap_mrc = mrcfile.open(emmap_path)
    apix = emmap_mrc.voxel_size.tolist()[0]
    origin = emmap_mrc.header.origin
    emmap = emmap_mrc.data
    
    
    # Takes the average of pixelsize in three axes for calculating filter cutoff
    filtered_emmap = low_pass_filter(emmap, dmin, apix)
    
    if output_filename is None:
        output_filename = emmap_path[:-4]+"_filtered_"+str(dmin)+"_A.mrc"
    
    save_as_mrc(filtered_emmap,output_filename=output_filename,apix=apix )
    
    return output_filename

def window3D(w):
    # Convert a 1D filtering kernel to 3D
    # eg, window3D(numpy.hanning(5))
    
    L=w.shape[0]
    m1=np.outer(np.ravel(w), np.ravel(w))
    win1=np.tile(m1,np.hstack([L,1,1]))
    m2=np.outer(np.ravel(w),np.ones([1,L]))
    win2=np.tile(m2,np.hstack([L,1,1]))
    win2=np.transpose(win2,np.hstack([1,2,0]))
    win=np.multiply(win1,win2)
    return win



def get_cosine_mask(mask,length_cosine_mask_1d):
    from scipy import signal
    cosine_window_1d = signal.cosine(length_cosine_mask_1d)
    cosine_window_3d = window3D(cosine_window_1d)
    cosine_mask = signal.convolve(mask,cosine_window_3d,mode='same')
    cosine_mask = cosine_mask/cosine_mask.max()
    return cosine_mask

def get_spherical_mask(mask_shape, radius_index):
    n = mask_shape[0]
    z,y,x = np.ogrid[-n//2:n//2,-n//2:n//2,-n//2:n//2]
    mask = (x**2+y**2+z**2 <= radius_index**2).astype(int)
    return mask

   
def resolution_cutoff(emmap_vol,cutoff_resolution,apix):
    nyquist_freq = 1/(apix*2)
    cutoff_freq = 1/(cutoff_resolution)
    n = emmap_vol.shape[0]
    index_of_cutoff_freq = round((n//2) * cutoff_freq/nyquist_freq)
    
    mask = get_spherical_mask(emmap_vol.shape,index_of_cutoff_freq)
    cosine_mask = get_cosine_mask(mask,length_cosine_mask_1d=n//10)
    cosine_mask = cosine_mask / cosine_mask.max()
    
    emmap_fourier_space = np.fft.fftn(emmap_vol)
    emmap_filtered = emmap_fourier_space * cosine_mask
    emmap_filtered_real_space = abs(np.fft.ifftn(emmap_filtered))
    
    return emmap_filtered_real_space

    

def butter_lowpass(cutoff, nyq_freq, order=4):
    from scipy import signal

    normal_cutoff = float(cutoff) / nyq_freq
    b, a = signal.butter(order, normal_cutoff, btype='lowpass')
    return b, a

def butter_lowpass_filter(data, cutoff_freq, nyq_freq, order=4):
    # Source: https://github.com/guillaume-chevalier/filtering-stft-and-laplace-transform
    from scipy import signal

    b, a = butter_lowpass(cutoff_freq, nyq_freq, order=order)
    y = signal.filtfilt(b, a, data)
    return y

       
def moving_average(array,window=5):
    return np.convolve(array,np.ones(window), 'same')/window