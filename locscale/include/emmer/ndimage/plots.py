#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 19 13:53:12 2022

@author: alok
"""

import numpy as np


def plot_radial_profile(freq,list_of_profiles, squared_amplitudes=False, legends=None, normalise=False, font=28,showlegend=True, showPoints=True, alpha=0.05, variation=None, yticks=None, logScale=True, ylims=None, xlims=None, crop_freq=None):
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
        
    i = 0
    colors = cm.rainbow(np.linspace(0,1,len(list_of_profiles)))
    

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
                if showPoints:
                    ax1.plot(freq**2,np.log(profile),c=colors[i], linewidth=1, marker="o")
                else:
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
                if showPoints:
                    ax1.plot(freq,profile,c=colors[i], linewidth=1, marker="o")
                else:
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
    if yticks is not None:
        plt.yticks(yticks)
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
    from locscale.include.emmer.ndimage.profile_tools import compute_radial_profile
    
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
    from locscale.include.emmer.ndimage.profile_tools import compute_radial_profile, frequency_array, plot_radial_profile, estimate_bfactor_through_pwlf
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