import numpy as np
#from fsc_attributes_calc import *
#from spline_fit_fsc import *
from fft_util import pyfftw_flag

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

def plot_fscs(dict_points,outfile,xlabel=None,ylabel=None,map_apix=1.0,
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
