#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec  3 15:58:13 2021

@author: alok
"""

## Plot tools 
import numpy as np
import seaborn as sns


    
def r2(y_fit, y_data):
    y_mean = y_data.mean()
    residual_squares = (y_data-y_fit)**2
    variance = (y_data-y_mean)**2
    
    residual_sum_of_squares = residual_squares.sum()
    sum_of_variance = variance.sum()
    
    r_squared = 1 - residual_sum_of_squares/sum_of_variance
    
    return r_squared

def crop_data_to_map(input_data_map, mask, mask_threshold, skip_zeros=True):
    from locscale.include.emmer.ndimage.map_utils import parse_input
    
    input_data_map = parse_input(input_data_map)
    mask = parse_input(mask)
    
    binarised_mask = (mask>=mask_threshold).astype(np.int_)
    flattend_array = (binarised_mask * input_data_map).flatten()
    
    nonzero_array = flattend_array[flattend_array>0]
    
    return nonzero_array

def plot_correlations(x_array, y_array, scatter=False, figsize=(14,8),font="Helvetica",fontscale=4,hue=None, x_label=None, y_label=None, title_text=None, output_folder=None, filename=None, find_correlation=True, alpha=0.3):
    
    import matplotlib as mpl
    import seaborn as sns
    import os
    import matplotlib.pyplot as plt
    from scipy import stats
    import pandas as pd

    import matplotlib as mpl
    mpl.rcParams['pdf.fonttype'] = 42
    sns.set(rc={'figure.figsize':figsize})
    sns.set_theme(context="paper", font=font, font_scale=fontscale)
    sns.set_style("white")
    
    fig = plt.figure(1)
    
    if x_label is None:
        x_label = "x"
    
    if y_label is None:
        y_label = "y"
    
    
    if find_correlation:
        data = pd.DataFrame(data=[x_array,y_array], index=[x_label,y_label]).T
      
        def annotate(data, **kws):
            r, p = stats.pearsonr(data[x_label], data[y_label])
            ax = plt.gca()
            ax.text(.05, .8, 'R$^2$={:.2f}'.format(r),
                    transform=ax.transAxes)
        g = sns.lmplot(data=data, x=x_label, y=y_label, scatter=scatter)
        #g.map_dataframe(annotate)
        plt.tight_layout()
        #plt.show()
    
    
        
        
    
    
    plt.xlabel(x_label)
    plt.ylabel(y_label)

    plt.tight_layout()
    if filename is not None:
        if output_folder is None:
            figure_output_folder = os.getcwd()
        else:
            figure_output_folder = output_folder
    
    
        output_filename = os.path.join(figure_output_folder, filename)
        plt.savefig(output_filename, dpi=600, bbox_inches="tight",transparency=True)
    else:
        return fig
        
def plot_correlations_multiple(xy_tuple, scatter=False, hue=None, x_label=None, y_label=None, ylims=None, title_text=None, output_folder=None, filename=None, find_correlation=True, alpha=0.3, ci=95):
    import seaborn as sns
    import os
    import matplotlib.pyplot as plt
    from scipy import stats
    import pandas as pd
    import matplotlib as mpl
    mpl.rcParams['pdf.fonttype'] = 42
    
    fig = plt.figure(figsize=(18,12))
    sns.set_theme(context="paper", font="Helvetica", font_scale=1.5)
    sns.set_style("white")
    kwargs = dict(linestyle="--", marker="o", linewidth=3, markersize=12, alpha=alpha)
    

    if x_label is None:
        x_label = "x"
    
    if y_label is None:
        y_label = "y"
    
    data_dictionary={}
    all_stacks = []
    for xy in xy_tuple:
        x_array, y_array, category_label = xy
        category_array = np.repeat(category_label, len(x_array))
        stack = np.vstack((x_array,y_array,category_array)).T
        all_stacks.append(stack)
    
    stack_arrays = np.concatenate(tuple([x for x in all_stacks]))
    
    data = pd.DataFrame(data=stack_arrays, columns=[x_label,y_label, "Category"])
    data[x_label] = data[x_label].astype(np.float32)
    data[y_label] = data[y_label].astype(np.float32)
    data["Category"] = data["Category"].astype(str)

        

    def annotate(data, **kws):
        r, p = stats.pearsonr(data[x_label], data[y_label])
        ax = plt.gca(figsize=(16,8))
        ax.text('R$^2$={:.2f}'.format(r),
                        transform=ax.transAxes)
    g = sns.lmplot(data=data, x=x_label, y=y_label, scatter=scatter, hue="Category", ci=ci, legend=False)
    
    plt.legend(loc="lower right")   

    plt.xlabel(x_label)
    plt.ylabel(y_label)
    
    if ylims is not None:
        plt.ylim(ylims)
    
        
    plt.tight_layout()
    if filename is not None:
        if output_folder is None:
            figure_output_folder = os.getcwd()
        else:
            figure_output_folder = output_folder
    
    
        output_filename = os.path.join(figure_output_folder, filename)
        plt.savefig(output_filename, dpi=600,bbox_inches='tight')
    else:
        return fig
    

def plot_linear_regression(data_input, x_col, y_col, x_label=None, y_label=None, title_text=None):
    import matplotlib.pyplot as plt
    import pandas as pd
    
    def linear(x,a,b):
        return a * x + b
    from matplotlib.offsetbox import AnchoredText
    f, ax = plt.subplots(1,1)

            
    data_unsort = data_input.copy()
    data=data_unsort.sort_values(by=x_col)
    x_data = data[x_col]
    y_data = data[y_col]
    
    y_fit = x_data ## When assuming y=x as ideal equation
    
    r_squared = data[x_col].corr(data[y_col], method="spearman")
    
    
    
    
    ax.plot(x_data, y_data,'bo')
    ax.plot(x_data, x_data, 'r-')
    equation = "y = x \nCorrelation = {}".format(round(r_squared,2))
    legend_text = equation
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
    
def plot(plot_properties):
    import seaborn as sns
    import matplotlib.pyplot as plt
    
    sns.set_theme(context="paper", font="Helvetica", font_scale=1.5)
    sns.set_style("white")
    kwargs = dict(linewidth=3)    
    

def plot_radial_profile_seaborn(freq, list_of_profiles, font=16, ylims=None, crop_first=10, crop_end=1, legends=None):
    import seaborn as sns
    import matplotlib.pyplot as plt
    
    import matplotlib as mpl
    mpl.rcParams['pdf.fonttype'] = 42
    
    freq = freq[crop_first:-crop_end]
    
    sns.set_theme(context="paper", font="Helvetica", font_scale=1.5)
    sns.set_style("white")
    kwargs = dict(linewidth=3)

    profile_list = np.array(list_of_profiles)
    average_profile = np.einsum("ij->j", profile_list) / len(profile_list)

    variation = []
    for col_index in range(profile_list.shape[1]):
        col_extract = profile_list[:,col_index]
        variation.append(col_extract.std())

    variation = np.array(variation)
        
    y_max = average_profile + variation
    y_min = average_profile - variation
    
    fig, ax = plt.subplots()
    ax = sns.lineplot(x=freq, y=average_profile[crop_first:-crop_end], **kwargs)
    ax.fill_between(freq, y_min[crop_first:-crop_end], y_max[crop_first:-crop_end], alpha=0.3)
    ax.set_xlabel('Spatial Frequency $1/d [\AA^{-1}]$',fontsize=font)
    ax.set_ylabel('$\mid F \mid $',fontsize=font)
    if legends is not None:
        ax.legend(legends)
    ax2 = ax.twiny()
    ax2.set_xticks(ax.get_xticks())
    ax2.set_xbound(ax.get_xbound())
    ax2.set_xticklabels([round(1/x,1) for x in ax.get_xticks()])
    ax2.set_xlabel('$d [\AA]$',fontsize=font)
    if ylims is not None:
        plt.ylim(ylims)
    plt.tight_layout()
    plt.show()
    
    return fig

def pretty_lineplot_XY(xdata, ydata, xlabel, ylabel, filename=None,figsize=(14,8), marker="o", markersize=12,fontscale=2.5,font="Helvetica",linewidth=2,legends=None):
    import matplotlib.pyplot as plt
    from matplotlib.pyplot import cm
    import seaborn as sns    
    from matplotlib.pyplot import cm
    import matplotlib as mpl
    ## Function not generic
    mpl.rcParams['pdf.fonttype'] = 42
    
    sns.set(rc={'figure.figsize':figsize})
    sns.set_theme(context="paper", font=font, font_scale=fontscale)
    sns.set_style("white")

    sns.lineplot(x=xdata,y=ydata,linewidth=linewidth,marker=marker,markersize=markersize)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)

    if legends is not None:        
        plt.legend(legends)
    plt.tight_layout()
    plt.savefig(filename, dpi=600, bbox_inches="tight", transparency=True)

def pretty_lineplot_multiple_fsc_curves(fsc_arrays_perturb, two_xaxis=True, filename=None,figsize=(14,8),fontscale=2.5,font="Helvetica",linewidth=2,legends=None):
    import matplotlib.pyplot as plt
    from matplotlib.pyplot import cm
    import seaborn as sns    
    from matplotlib.pyplot import cm
    import matplotlib as mpl
    ## Function not generic
    mpl.rcParams['pdf.fonttype'] = 42
    fig = plt.figure()
    sns.set(rc={'figure.figsize':figsize})
    sns.set_theme(context="paper", font=font, font_scale=fontscale)
    sns.set_style("white")
    fsc_filename = filename
    colors_rainbow = cm.rainbow(np.linspace(0,1,len(fsc_arrays_perturb.keys())))
    
    if two_xaxis:
        # print(';)
        fig = plt.figure()
        ax1 = fig.add_subplot(111)
        ax1.grid(False)
        
        
        for i,rmsd in enumerate(fsc_arrays_perturb.keys()):
            sns.lineplot(x=fsc_arrays_perturb[rmsd][0],y=fsc_arrays_perturb[rmsd][1], linewidth=linewidth, color=colors_rainbow[i], ax=ax1)
            ax1.set_xlabel(r" Spatial Frequency, $d^{-1}(\AA^{-1}$)")
            ax1.set_ylabel("FSC")
        
        if legends is not None:        
            ax1.legend(legends)
        
        ax2 = ax1.twiny()
        ax2.set_xticks(ax1.get_xticks())
        ax2.set_xbound(ax1.get_xbound())        
        ax2.set_xticklabels([round(1/x,1) for x in ax1.get_xticks()])            
        ax2.set_xlabel(r'Resolution, $d (\AA)$')
        
        if legends is not None:   
            print("Legends print")
        plt.legend(legends)
    else:
        for i,rmsd in enumerate(fsc_arrays_perturb.keys()):
            sns.lineplot(x=fsc_arrays_perturb[rmsd][0],y=fsc_arrays_perturb[rmsd][1], linewidth=linewidth, color=colors_rainbow[i])
            plt.xlabel(r" Spatial Frequency, $d^{-1}(\AA^{-1}$)")
            plt.ylabel("FSC")
    

    #plt.tight_layout()
    fig.savefig(fsc_filename, dpi=600, bbox_inches="tight", transparency=True)

def pretty_violinplots(list_of_series, xticks, ylabel,xlabel=None, figsize=(14,8),fontscale=3,font="Helvetica",linewidth=2, filename=None):
    import seaborn as sns
    import matplotlib.pyplot as plt
    from matplotlib.pyplot import cm
    import matplotlib as mpl
    ## Function not generic
    mpl.rcParams['pdf.fonttype'] = 42
    fig = plt.figure()
    sns.set(rc={'figure.figsize':figsize})
    sns.set_theme(context="paper", font=font, font_scale=fontscale)
    sns.set_style("white")
    
    ax = sns.violinplot(data=list_of_series, scale_hue=False)
    ax.set_xticklabels(xticks)
    ax.set_ylabel(ylabel)
    if xlabel is not None:
        ax.set_xlabel(xlabel)
    fig.tight_layout()
    
    if filename is not None:
        plt.savefig(filename, dpi=600, bbox_inches="tight", transparency=True)
        

def pretty_boxplots(list_of_series, xticks, ylabel,xlabel=None, figsize=(14,8),fontscale=3,font="Helvetica",linewidth=2, filename=None):
    import seaborn as sns
    import matplotlib.pyplot as plt
    from matplotlib.pyplot import cm
    import matplotlib as mpl
    ## Function not generic
    mpl.rcParams['pdf.fonttype'] = 42
    fig,ax = plt.subplots()
    sns.set(rc={'figure.figsize':figsize})
    sns.set_theme(context="paper", font=font, font_scale=fontscale)
    sns.set_style("white")
    
    ax.boxplot(list_of_series)
    ax.set_xticklabels(xticks)
    ax.set_ylabel(ylabel)
    if xlabel is not None:
        ax.set_xlabel(xlabel)
    fig.tight_layout()
    
    if filename is not None:
        plt.savefig(filename, dpi=600, bbox_inches="tight", transparency=True)
    
    
def pretty_plot_radial_profile(freq,list_of_profiles_native,normalise=True, discrete=True, legends=None,figsize=(14,8),fontsize=28, linewidth=1, marker="o", font="Helvetica",fontscale=2.5, showlegend=True, showPoints=False, alpha=0.05, variation=None, yticks=None, logScale=True, ylims=None, xlims=None, crop_freq=None):
    import matplotlib.pyplot as plt
    from matplotlib.pyplot import cm
    from locscale.include.emmer.ndimage.profile_tools import crop_profile_between_frequency
    import seaborn as sns

    import matplotlib as mpl
    mpl.rcParams['pdf.fonttype'] = 42
    
    sns.set(rc={'figure.figsize':figsize})
    sns.set_theme(context="paper", font=font, font_scale=fontscale)
    sns.set_style("white")
    
    if normalise:
        list_of_profiles = []
        for profile in list_of_profiles_native:
            normalised_profile = profile/profile.max()
            list_of_profiles.append(normalised_profile)
    else:
        list_of_profiles = list_of_profiles_native
        

    i = 0
    colors = cm.rainbow(np.linspace(0,1,len(list_of_profiles)))
    
    if discrete:
        fig = plt.figure()
        ax1 = fig.add_subplot(111)
        ax1.grid(False)
        ax2 = ax1.twiny()

        if logScale:
            for profile in list_of_profiles:
                if crop_freq is not None:
                    frequency, profile = crop_profile_between_frequency(freq, profile, crop_freq[0], crop_freq[1])
                if showPoints:
                    ax1.plot(frequency**2,np.log(profile),c=colors[i], linewidth=linewidth, marker=marker)
                else:
                    ax1.plot(frequency**2,np.log(profile),c=colors[i], linewidth=linewidth)
                i += 1
            
            ax2.set_xticks(ax1.get_xticks())
            ax2.set_xbound(ax1.get_xbound())
            ax2.set_xticklabels([round(1/np.sqrt(x),1) for x in ax1.get_xticks()])
            if showlegend:
                ax1.legend(legends)
            ax1.set_xlabel(r'Spatial Frequency, $d^{-2} (\AA^{-2})$')
            ax1.set_ylabel(r'$ln  \langle \mid F \mid \rangle $ ')
            ax2.set_xlabel(r'Resolution, $d (\AA)$')
        else:
            for profile in list_of_profiles:
                if crop_freq is not None:
                    frequency, profile = crop_profile_between_frequency(freq, profile, crop_freq[0], crop_freq[1])
                if showPoints:
                    ax1.plot(frequency,profile,c=colors[i], linewidth=linewidth, marker="o")
                else:
                    ax1.plot(frequency,profile,c=colors[i], linewidth=linewidth)
                i += 1
            
            
            if showlegend:
                ax1.legend(legends)
        
                
            ax2.set_xticks(ax1.get_xticks())
            ax2.set_xbound(ax1.get_xbound())
            ax2.set_xticklabels([round(1/x,1) for x in ax1.get_xticks()])
            
    
            ax1.set_xlabel(r'Spatial Frequency, $d^{-1} [\AA^{-1}]$')
            ax1.set_ylabel(r'Normalised $ \langle F \rangle $')
            ax2.set_xlabel(r'Resolution, $d (\AA)$')

    else:
        
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
                frequency, average_profile = crop_profile_between_frequency(freq, average_profile, crop_freq[0], crop_freq[1])
                frequency, y_max = crop_profile_between_frequency(freq, y_max, crop_freq[0], crop_freq[1])
                frequency, y_min = crop_profile_between_frequency(freq, y_min, crop_freq[0], crop_freq[1])
            
            ax1.plot(frequency**2, np.log(average_profile), 'k',alpha=1)
            ax1.fill_between(frequency**2,np.log(y_max), np.log(y_min), color="grey", alpha=0.5)
            if showlegend:
                ax1.legend(["N={}".format(len(profile_list))])
        
            
            ax2.set_xticks(ax1.get_xticks())
            ax2.set_xbound(ax1.get_xbound())
            ax2.set_xticklabels([round(1/np.sqrt(x),1) for x in ax1.get_xticks()])
            
    
            ax1.set_xlabel(r'Spatial Frequency, $d^{-2} [\AA^{-2}]$')
            ax1.set_ylabel(r'$ln  \langle \mid F \mid \rangle $ ')
            ax2.set_xlabel(r'Resolution, $d (\AA)$')
        else:
            if crop_freq is not None:
                frequency, average_profile = crop_profile_between_frequency(freq, average_profile, crop_freq[0], crop_freq[1])
                frequency, y_max = crop_profile_between_frequency(freq, y_max, crop_freq[0], crop_freq[1])
                frequency, y_min = crop_profile_between_frequency(freq, y_min, crop_freq[0], crop_freq[1])
            ax1.plot(frequency, average_profile, 'k',alpha=1)
            ax1.fill_between(frequency,y_max, y_min,color="grey", alpha=0.5)
            
            if showlegend:
                ax1.legend(["N={}".format(len(profile_list))])
        
                
            ax2.set_xticks(ax1.get_xticks())
            ax2.set_xbound(ax1.get_xbound())
            ax2.set_xticklabels([round(1/x,1) for x in ax1.get_xticks()])
            
    
            ax1.set_xlabel(r'Spatial Frequency, $d^{-1} [\AA^{-1}]$')
            ax1.set_ylabel(r'Normalised $ \langle F \rangle $')
            ax2.set_xlabel(r'Resolution, $d (\AA)$')

    if ylims is not None:
        plt.ylim(ylims)
    if yticks is not None:
        plt.yticks(yticks)
    if xlims is not None:
        plt.xlim(xlims)
    
    
    plt.tight_layout()
    return fig
    