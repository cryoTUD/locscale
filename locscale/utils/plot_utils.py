#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec  3 15:58:13 2021

@author: alok
"""

## Plot tools 

def r2(y_fit, y_data):
    y_mean = y_data.mean()
    residual_squares = (y_data-y_fit)**2
    variance = (y_data-y_mean)**2
    
    residual_sum_of_squares = residual_squares.sum()
    sum_of_variance = variance.sum()
    
    r_squared = 1 - residual_sum_of_squares/sum_of_variance
    
    return r_squared

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
    
    
    
    