import sys
import fft_calc
from spline_fit_fsc import *
from scipy import interpolate, optimize
from scipy.signal import argrelextrema
import numpy as np
import random
import warnings

def get_resolution_from_fsc_curve(listfreq,listfsc,map_apix,
                                  fsc_cutoff=0.5,minRes=20.0,
                                  maxRes=1.5,returnsols=True,
                                  list_weights=None):
    '''
    Get resolution at FSC 0.5 from the freq v FSC plot
    
    '''
    solutions = []
    #f = interpolate.interp1d(list1, list2, kind='linear',fill_value='extrapolate')
    #try: fsc05 = f(0.5)
    #except ValueError: return 0.0
    #get roots of the fsc curve at 0.5 using spline fit
    solutions, listspl = get_roots_from_fsccurve(
                                                listfreq,listfsc,
                                                fsc_cutoff=fsc_cutoff,
                                                list_weights=list_weights)
    #polynomial doesnt cross 0.5
    if len(solutions) == 0:
        if listfsc[-1] > fsc_cutoff and listfsc[0] > fsc_cutoff:
            #find roots based on local minima
            fsc05 = get_fsc_based_on_local_minima(listfreq,listfsc,
                                                    map_apix,
                                                    maxRes=maxRes,
                                                    minRes=minRes)
            solutions = [fsc05]
        else:
            solutions = [map_apix/minRes]
            fsc05 = map_apix/minRes
        
    else:
        #check if sol is in fall off of curve
        fsc05 = select_fsc_from_roots(solutions,listfreq,listspl,
                                      map_apix,fsc_cutoff=fsc_cutoff,
                                      maxRes=maxRes,minRes=minRes)
        solutions = [fsc05]
    try: resolution = map_apix/fsc05
    except UnboundLocalError: resolution = minRes
    #if returnsols: return resolution, solutions
    return resolution, solutions

def select_fsc_from_roots(solutions, listfreq,listfsc,map_apix,
                          fsc_cutoff=0.5,minRes=20.0,maxRes=1.5):
    '''
    Select a FSC0.5 from list of solutions to the fsc curve
    
    '''
    sorted_sol = np.sort(solutions)
    list_selfsc05 = []
    list_slopes = []
    list_fsc05 = []
    for i in range(1,len(solutions)+1):
        #set res as minRes if no solutions are better
        if map_apix/sorted_sol[-i] > minRes:
            fsc05 = map_apix/minRes
#             if len(list_selfsc05) == 0: return fsc05
            if len(list_fsc05) == 0: 
                list_fsc05.append(fsc05)
            break
        #set res as maxRes if no solutions are better
        elif map_apix/sorted_sol[-i] < maxRes:
            if len(solutions) > 1: 
                try:
                    #continue to next solution if one exists within the range
                    if map_apix/sorted_sol[-(i+1)] < minRes:
                        continue
                except IndexError: return map_apix/maxRes
            else: return map_apix/maxRes
        fsc05 = sorted_sol[-i]
        try:
            #skip possible (drop and) rise around FSC 0.5
            if map_apix/sorted_sol[-(i+1)] - map_apix/sorted_sol[-i] < 3.0:
                continue
        except IndexError:
            pass
        list_fsc05.append(fsc05)

#         index = np.searchsorted(listfreq,fsc05,side='right')
#         try:
#             #prev point
#             if listfsc[index-1]>fsc_cutoff:
#                 flag_prevpoint=0
#                 #two points before
#                 if index-2 >= 0:
#                     if listfsc[index-2]>fsc_cutoff: 
#                         slope = (listfsc[index-2] - fsc_cutoff)/(fsc05-listfreq[index-2])
#                         list_selfsc05.append(fsc05)
#                         list_slopes.append(slope)
#                         flag_prevpoint = 1
#                 if index+1 < len(listfsc):
#                     if listfsc[index+1]< fsc_cutoff: 
#                         slope_nextpoint = (listfsc[index-1] - listfsc[index+1])/  \
#                                 (listfreq[index+1] - listfreq[index-1])
#                         if flag_prevpoint == 1:
#                             list_slopes[-1] = max(slope,slope_nextpoint)
#                         else:
#                             list_selfsc05.append(fsc05)
#                             list_slopes.append(slope_nextpoint)
#                             
#                 else: break
# #             elif listfsc[index+1]< fsc_cutoff:
# #                 if index+2 < len(listfsc):
# #                     if listfsc[index+2]< fsc_cutoff: break
#         except IndexError: break
#     if len(list_slopes) > 0:
#         fsc05 = list_selfsc05[list_slopes.index(max(list_slopes))]
    #return solution close to maximum resolution
    return list_fsc05[-1]

def select_min_and_max_freq(sols,list_freq,n=100,minperbin=5,
                            minFreq=0.0,maxFreq=0.5):
    '''
    Identify min and max freq bounds for fsc 0.5 from the sample set
    '''
    #calculate histogram of fsc05 freqs
    freq,bins = np.histogram(sols,list_freq)
    freq_step = list_freq[1]-list_freq[0]
    sum_freq = 0
    dict_freq = {}
    not_minperbin = 0
    for i in range(1,len(freq)+1):
        #from high to low freq
        f = freq[-i]
        if f > minperbin:
            #prev bin had less than minperbin?
            if sum_freq == 0: 
                #set end of bin
                binend = bins[-i]
                #reset flag
                not_minperbin = 0
            
            binstart = bins[-i-1]
            sum_freq += f
        else:
            #if prev bins had more than minperbin and
            #last bin had less than minperbin (to skip any lone bins between 
            #bins with solutions)
            if sum_freq > 0 and not_minperbin == 1: 
                #save prev bin bounds
                fraction = round(float(sum_freq)/n,2)
                if fraction > 0.50:
                    minfreq_sel, maxfreq_sel = \
                            max(binstart-0.05*freq_step,minFreq),\
                            min(binend+0.05*freq_step,maxFreq)
                    return minfreq_sel, maxfreq_sel
                #reset sum_freq
                sum_freq = 0
            not_minperbin = not_minperbin+1
    if sum_freq > 0: 
        fraction = round(float(sum_freq)/n,2)
        if fraction > 0.50:
            minfreq_sel, maxfreq_sel = \
                    max(binstart-0.05*freq_step,minFreq),\
                    min(binend+0.05*freq_step,maxFreq)
            return minfreq_sel, maxfreq_sel
        
    print 'Potential artifacts in the local FSC curves'
    warnings.warn("Check for map artifacts/sytematic noise (e.g. tight mask)")
    return None, None

def stack_1d_arrays(inparray,onedarray):
    return np.vstack((inparray,onedarray))

def extend_minfreq_by_deviation(dict_points,minRes,maxRes,apix):
    """
    Extend the minimum frequency bound based on FSC deviation in each shell
    """
    minfreq = apix/minRes
    i = 0
    for k in dict_points:
        if i == 0:
            listfreq = dict_points[k][0]
            stacked_array = dict_points[k][1]
        else: stacked_array = stack_1d_arrays(stacked_array,dict_points[k][1])
        i += 1
    mean_array = np.array([np.mean(stacked_array[:,i]) \
                    for i in range(len(stacked_array[0]))])
    #indices of local minima
    local_min = argrelextrema(mean_array, np.less)[0]
    minFreqSel = None
    maxFreqSel = None
    for indi in local_min:
        if indi > 0:
            if mean_array[indi] < 0.8 and listfreq[indi] < minfreq:
                try: 
                    #save the first local minima
                    if minFreqSel is None: 
                        minFreqSel = listfreq[indi+1]
                        continue
                except IndexError: pass
            elif abs(apix/listfreq[indi] - maxRes) < 5.0:
                #save the last minima close to maxRes
                maxFreqSel = listfreq[indi]
    if not minFreqSel is None:
        return minFreqSel, maxFreqSel
    mad = np.array([np.median(np.absolute(
                    stacked_array[:,i]-np.median(stacked_array[:,i]))) \
                    for i in range(len(stacked_array[0]))])
    
    median_z = mad/np.median(mad)
    
    minfreq_ind = np.searchsorted(listfreq,minfreq)
    i = 0
    
    #first index not included, minimum index is 2 (i+1)
    for i in np.arange(max(0,minfreq_ind-1),-1,-1):
        if median_z[i] > max(1.5,median_z[minfreq_ind]) : break
    if i == 0: minfreq = listfreq[0]
    elif i > 0: minfreq = listfreq[i+1]
    return minfreq, None
#     mad = (mad - min(mad)) / (max(mad) - min(mad))
#     plt.plot(x_array,mad,linewidth=3.0,marker='*',linestyle='--',
#              color='black')

def get_medians_fsc(dict_points):
    '''
    get medians based on fsc distribution at each shell
    dict_points is a dict with freq and fsc lists from a sample set
    '''
    i = 0
    for k in dict_points:
        if i == 0:
            #listfreq = dict_points[k][0]
            stacked_array = dict_points[k][1]
        else: stacked_array = stack_1d_arrays(stacked_array,dict_points[k][1])
        i += 1
    list_medians = [np.median(stacked_array[:,i]) 
                    for i in range(len(stacked_array[0]))]
    mad = [np.median(np.absolute(
                    stacked_array[:,i]-np.median(stacked_array[:,i]))) \
                    for i in range(len(stacked_array[0]))]
    return list_medians, mad

def calculate_weights_from_median_deviation(listfsc,list_medians,
                                            mad=None):
    '''
    calculate weights [0 - 1.] for fsc values in each shell
    '''
    if list_medians is None: return np.ones(len(listfsc))
    #deviation from median
#     dev = np.array([np.absolute(listfsc[i]-list_medians[i]) 
#                     for i in range(len(listfsc))])
#     dev1 = np.array(mad)/dev
    dev = np.ones(len(listfsc))
    #set low weights for last 3 values
    dev[-3:] = 0.1
    for i in range(1,len(listfsc)):
        #include all high fsc values?
        if listfsc[i] > 0.75: continue
        #compare median trend with trend in current fsc
        if listfsc[i-1]-listfsc[i] == 0.0: continue
        prevdiff_ratio = (list_medians[i-1]-list_medians[i])/ \
                            (listfsc[i-1]-listfsc[i])
        #opposite to the median slope?
        if prevdiff_ratio < -1.:  
            dev[i] = 0.0
            #reduce the next as well
            if i+1 < len(dev): 
                dev[i+1] = 0.3
        #significantly higher compared to median trend
        elif prevdiff_ratio > 2.0: 
            dev[i] = 0.0
    #high freq > 0.5?
    for i in xrange(3):
        if listfsc[-i] > 0.5:
            dev[-i] = 0.0
    for i in range(1,len(listfsc)-1):
        if listfsc[i] > 0.75: continue
        #isolated maxima
        if dev[i] == 1. and not dev[i-1] == 1. and not dev[i+1] == 1.:
            dev[i] = (dev[i-1]+dev[i+1])/2.
    #print dev
    return dev

def calc_area_under_fsc_curve(listfreq,listfsc,
                              map_apix,
                              maxRes=1.5,minRes=20.0,
                              list_weights=None):
    maxfreq = map_apix/maxRes
    minfreq = map_apix/minRes
    minfreq_index = np.searchsorted(listfreq,minfreq)
    maxfreq_index = np.searchsorted(listfreq,maxfreq,side='right')
    maxfreq_index = min(maxfreq_index,len(listfreq))
    #temporarily raise warnings as errors
    with warnings.catch_warnings():
        warnings.filterwarnings('error')
#     try:
#         listfsc_spline = get_spline_fit_fsc(listfreq,listfsc,
#                                         list_weights=list_weights)
#     except RuntimeWarning:
#         listfsc_spline = listfsc
#     listfsc_spline = listfsc
    auc = np.trapz(listfsc[minfreq_index:maxfreq_index],
             listfreq[minfreq_index:maxfreq_index])
    avg_auc = auc/(maxfreq_index-minfreq_index)
    #if avg_auc > 0.026: plot_spline(listfreq, listfsc, listfsc_spline)
    return avg_auc

def calculate_fscavg(listfreq,listfsc,listnsf,
                              map_apix,
                              maxRes=None,minRes=None,
                              list_weights=None):
    if maxRes is None:
        maxfreq_index = len(listfreq)
    else:
        maxfreq = map_apix/maxRes
        print maxfreq
        maxfreq_index = np.searchsorted(listfreq,maxfreq,side='right')
    if minRes is None:
        minfreq_index = 0
    else:
        minfreq = map_apix/minRes
        minfreq_index = np.searchsorted(listfreq,minfreq)
        
    maxfreq_index = min(maxfreq_index,len(listfreq))
    weighted_fsc_sum = 0.0
    sum_nsf = 0.
    for ind in xrange(minfreq_index,maxfreq_index):
        print listfreq[ind], listfsc[ind],listnsf[ind]
        weighted_fsc_sum += listfsc[ind]*listnsf[ind]
        sum_nsf += listnsf[ind] 
    FSCavg = weighted_fsc_sum/sum_nsf
    return FSCavg
    
def calc_z_scores_median(list_vals):
    array_vals = np.array(list_vals)
    diff_median = array_vals-np.median(list_vals)
    mad = np.median(np.absolute(array_vals-np.median(list_vals)))
    z_median = diff_median/mad
    return z_median
    
