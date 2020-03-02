import numpy as np
from scipy.signal import argrelextrema
from scipy import interpolate
import warnings

#get freq-FSC curve (cubic)
def get_spline_fit_fsc(listfreq,listfsc,
                       plot=True,
                       list_weights=None):
    #fit to cubic polynomial
    tck = fit_spline_xy_cubic(listfreq,listfsc,list_weights=list_weights)
    listspl = interpolate.splev(listfreq, tck)
    return listspl

def plot_spline(l1,l2,spl):
    import matplotlib.pyplot as plt
    plt.plot(l1, l2, 'g', lw=2,ls='-')
    plt.plot(l1, spl, 'g', lw=3,ls='--')
    plt.show()

def fit_spline_xy_cubic_extrapolate(listfreq,listfsc,fsc_cutoff=0.5):
    f = interpolate.UnivariateSpline(listfreq, np.array(listfsc)-fsc_cutoff, 
                                     k=3,
                                     ext=0)
    
    return f.root()

#fit freq-FSC curve (cubic)
def fit_spline_xy_cubic(listfreq,listfsc,
                        list_weights=None):
    #fit to cubic polynomial
    if list_weights is None: 
        tck = interpolate.splrep(listfreq, listfsc,k=3,s=0.005)
    else:
        tck = interpolate.splrep(listfreq, listfsc,w=list_weights,
                                 k=3,s=0.005)
    return tck

def get_roots_from_fsccurve(listfreq,listfsc,fsc_cutoff=0.5,
                        plot=False,
                        list_weights=None):
    '''
    Get roots of the fsc spline fit at 0.5
    Return:
        Roots, list of spline values
    '''
    #temporarily raise warnings as errors
    with warnings.catch_warnings():
        warnings.filterwarnings('error')
    try:
        tck = fit_spline_xy_cubic(listfreq,listfsc,list_weights=list_weights)
    except RuntimeWarning: return [], []
    tck_mod = (tck[0], tck[1] - fsc_cutoff, tck[2])
    #roots at fsc 0.5
    solutions = interpolate.sproot(tck_mod)
    try: listspl = interpolate.splev(listfreq, tck)
    except: listspl = listfsc
    if plot: 
        plot_spline(listfreq,listfsc,listspl)
    #check for raise of the spline curve at the end
    if listspl[-1] > listfsc[-1]: listspl[-1] = listfsc[-1]
    return solutions, listspl

def select_fsc_from_local_minimas(locminima_indices,sorted_freqarray,
                                   sorted_fscarray,map_apix,minRes=20.0,
                                   maxRes=1.5, fsc_cutoff=0.5):
    '''
    From all local minimas, 
    select one local minima within the max and min resolutions
    '''
    #This part doesnt seem to be used?
    try:
        #assuming fsc has 'one' minima around optimal resolution
        #set fsc0.5 as freq between last min and next freq 
        #this helps when the fsc line doesnt fall below 0.5
        fsc05 = (sorted_freqarray[locminima_indices[-1]]+
             sorted_freqarray[locminima_indices[-1]+1])/2.
    except IndexError: 
        #set as last minima
        fsc05 = sorted_freqarray[locminima_indices[-1]]
        
    #scan from right to left for minimas
    for i in range(1,len(locminima_indices)+1):
        lastminima = locminima_indices[-i]
        #if resolution at minima is beyond minimun resolution expected
        if map_apix/sorted_freqarray[lastminima] > minRes: 
            #set fsc05 as minimum expected resolution
            fsc05 = map_apix/minRes
            return fsc05
        #if resolution at minima is better than maximum expected resolution
        elif map_apix/sorted_freqarray[lastminima] < maxRes:
            try:
                #continue to lower frequencies if there are minimas in the expected range
                if map_apix/locminima_indices[-(i+1)] < minRes:
                    continue
            except IndexError: 
                return map_apix/maxRes
            #return fsc05 at maximum expected resolution
            return map_apix/maxRes
        #alternately select a neighboring freq to calculate average for fsc05
        try:
            #try freq at the right
            neighfreq = sorted_freqarray[lastminima+1]
            #if fsc curve raises from a lower value on the left to a higher value on the right
            if lastminima-1 >= 0:
                if sorted_fscarray[lastminima-1] < \
                     sorted_fscarray[lastminima+1]:
                    fsc05 = sorted_freqarray[lastminima]
                else:
                    #select the neighbor on the right
                    neighfreq = sorted_freqarray[lastminima-1]
                    fsc05 = (sorted_freqarray[lastminima]+neighfreq)/2.
            else:
                fsc05 = sorted_freqarray[lastminima]
            break
        except IndexError: pass
    return fsc05
    
def get_fsc_based_on_local_minima(listfreq,listfsc,map_apix,minRes=20.0,
                                  maxRes=1.5, fsc_cutoff=0.5):
    '''
    Find all local minimas and select one within expected resolution range
    '''
    fscarray = np.array(listfsc)
    freqarray = np.array(listfreq)
    #use local minima as fsc 0.5
    locminima = argrelextrema(fscarray, np.less,
                              order=3)[0]
    if len(locminima) != 0:
        fsc05 = select_fsc_from_local_minimas(locminima,freqarray,
                                               fscarray,
                                               map_apix=map_apix,
                                               minRes=minRes,
                                               maxRes=maxRes)
    #select highest freq
    else: fsc05 = listfreq[-1]
    return fsc05
