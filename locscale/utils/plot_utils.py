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

def plot_correlations(x_array, y_array, scatter=False, hue=None, x_label=None, y_label=None, title_text=None, output_folder=None, filename=None, find_correlation=True, alpha=0.3):
    import seaborn as sns
    import os
    import matplotlib.pyplot as plt
    from scipy import stats
    import pandas as pd
    fig = plt.figure(1)
    sns.set_theme(context="paper", font="Helvetica", font_scale=1.5)
    sns.set_style("white")
    kwargs = dict(linestyle="--", marker="o", linewidth=3, markersize=12, alpha=alpha)


    
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
        plt.savefig(output_filename, dpi=600)
    else:
        return fig
        
def plot_correlations_multiple(xy_tuple, scatter=False, hue=None, x_label=None, y_label=None, title_text=None, output_folder=None, filename=None, find_correlation=True, alpha=0.3, ci=95):
    import seaborn as sns
    import os
    import matplotlib.pyplot as plt
    from scipy import stats
    import pandas as pd
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
    
def default_specification():    
    plot_specification = {}
    plot_specification['dpi'] = 600
    plot_specification['figsize'] = (8,8)
    plot_specification['global_font'] = 'Helvetica'
    

def plot_vector_graphics(plot_data, filepath, plot_specification):
    import matplotlib.pyplot as plt
    
    x_data = plot_data['x_data']
    y_data = plot_data['y_data']
    xlabel = plot_data['x_label']
    ylabel = plot_data['y_label']
    title = plot_data['title']
    
    
    
    