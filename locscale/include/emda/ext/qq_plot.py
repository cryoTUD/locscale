"""
Author: "Rangana Warshamanage, Garib N. Murshudov"
MRC Laboratory of Molecular Biology
    
This software is released under the
Mozilla Public License, version 2.0; see LICENSE.
"""

# Quantile-Quantile plots for distributions (data) analysis

from __future__ import absolute_import, division, print_function, unicode_literals
import fcodes_fast
import numpy as np
from emda import core, ext
from emda.config import *
from matplotlib import pyplot as plt


def pass_cif_sf(cif_name):
    import gemmi

    doc = gemmi.cif.read(cif_name)
    rblock = gemmi.as_refln_blocks(doc)[0]
    # h = rblock.make_array_int('index_h', -1000)
    # k = rblock.make_array_int('index_k', -1000)
    # l = rblock.make_array_int('index_l', -1000)
    fobs = rblock.make_array_float("F_meas_au")
    resol = rblock.make_1_d2_array().round(4)
    resol = np.sqrt(1.0 / resol)
    return fobs, resol


def make_xycols(set1, set2, resol1, resol2):
    lowres = np.max([np.max(resol1), np.max(resol2)])
    higres = np.min([np.min(resol1), np.min(resol2)])
    res_arr = fcodes_fast.setbin(10, lowres, higres, binmode)
    bin_idx1 = fcodes_fast.binarr_1d(res_arr, resol1, len(res_arr), len(resol1))
    bin_idx2 = fcodes_fast.binarr_1d(res_arr, resol2, len(res_arr), len(resol2))
    xq_lst = []
    yq_lst = []
    for ibin in range(len(res_arr)):
        print(ibin)
        x = set1 * (bin_idx1 == ibin)
        x = x[~np.isnan(x)]  # remove NaN values
        # x = x[x > 0.]
        # if x.size > 0: x = x/np.max(x)
        # normalize to variance of Fobs
        x_variance = np.average((x - np.average(x)) ** 2)
        if x.size > 0:
            x = x / np.sqrt(x_variance)
        y = set2 * (bin_idx2 == ibin)
        y = y[~np.isnan(y)]
        # y = y[y > 0.]
        # if y.size > 0: y = y/np.max(y)
        y_variance = np.average((y - np.average(y)) ** 2)
        if y.size > 0:
            y = y / np.sqrt(y_variance)
        x_q, y_q = qqplot_simple(x, y)
        xq_lst.append(x_q)
        yq_lst.append(y_q)
        print(x_variance, y_variance)
    return xq_lst, yq_lst


def using_3d_data(uc, f1, f2):
    # create bin_idx
    nx, ny, nz = f1.shape
    maxbin = np.amax(np.array([nx // 2, ny // 2, nz // 2]))
    nbin, _, bin_idx = fcodes_fast.resolution_grid(uc, debug_mode, maxbin, nx, ny, nz)
    I1 = np.real(f1 * f1.conj())
    I2 = np.real(f2 * f2.conj())
    xq_lst = []
    yq_lst = []
    for ibin in range(15, 20):
        x = (I1[bin_idx == ibin]).flatten()
        x = x[~np.isnan(x)]  # remove NaN values
        x_variance = np.average((x - np.average(x)) ** 2)
        # normalise to sqrt of variance
        if x.size > 0:
            x = x / np.sqrt(x_variance)
        y = (I2[bin_idx == ibin]).flatten()
        y = y[~np.isnan(y)]
        y_variance = np.average((y - np.average(y)) ** 2)
        if y.size > 0:
            y = y / np.sqrt(y_variance)
        xq_lst.append(x)
        yq_lst.append(y)
        print(x_variance, y_variance)
    return xq_lst, yq_lst


def threed_data(uc, f1, f2):
    assert f1.shape == f2.shape
    nx, ny, nz = f1.shape
    f1d_1, res1d_1 = fcodes_fast.conv3d_to_1d(f1, uc, 1, nx, ny, nz)
    f1d_2, res1d_2 = fcodes_fast.conv3d_to_1d(f2, uc, 1, nx, ny, nz)
    xq_lst, yq_lst = make_xycols(f1d_1, f1d_2, res1d_1, res1d_2)
    return xq_lst, yq_lst


def draw_qqplot(xq_lst, yq_lst):
    fig = plt.figure(figsize=(6, 6))
    ax1 = fig.add_subplot(111)
    ibin = -1
    xmax = 0
    ymax = 0
    for x_quantiles, y_quantiles in zip(xq_lst, yq_lst):
        ibin = ibin + 1
        if x_quantiles.size > 0 and xmax < np.max(x_quantiles):
            xmax = np.max(x_quantiles)
        if y_quantiles.size > 0 and ymax < np.max(y_quantiles):
            ymax = np.max(y_quantiles)
        ax1.scatter(
            x_quantiles,
            y_quantiles,
            s=10,
            marker="o",
            label="bin_" + str(ibin),
            alpha=0.3,
            edgecolors="none",
        )
    plt.title("Q-Q Plot of Fobs of 6pu5(1st) vs 6pu4(2nd)")
    plt.xlabel("Quantiles of 1st sample")
    plt.ylabel("Quantiles of 2nd sample")
    plt.plot((0, ymax), (0, ymax), "lightgray")
    plt.legend(loc=0)
    plt.savefig("qqplot.eps", format="eps", dpi=300)
    plt.show()
    plt.close()


def make_quantiles(x, y):
    import statsmodels.api as sm
    from statsmodels.graphics.gofplots import qqplot_2samples

    pp_x = sm.ProbPlot(x)
    pp_y = sm.ProbPlot(y)
    qqplot_2samples(pp_x, pp_y, line="45")
    plt.show()


def qqplot_simple(x, y, quantiles=None, interpolation="nearest"):
    """Some insights: 
    https://stats.stackexchange.com/questions/403652/two-sample-quantile-quantile-plot-in-python"""

    if quantiles is None:
        quantiles = min(len(x), len(y))
    # Compute quantiles of the two samples
    if isinstance(quantiles, int):
        quantiles = np.linspace(start=0, stop=1, num=int(quantiles))
    else:
        quantiles = np.atleast_1d(np.sort(quantiles))
    x_quantiles = np.quantile(x, quantiles, interpolation=interpolation)
    y_quantiles = np.quantile(y, quantiles, interpolation=interpolation)
    return x_quantiles, y_quantiles


def qqplot(input_map1, input_map2):
    if input_map1.endswith((".mrc", ".map")):
        uc, arr1, _ = core.iotools.read_map(input_map1)
        uc, arr2, _ = core.iotools.read_map(input_map2)
        f1 = np.fft.fftshift(np.fft.fftn(arr1))
        f2 = np.fft.fftshift(np.fft.fftn(arr2))
        # xq_lst, yq_lst = using_3d_data(uc, f1, f2)
        xq_lst, yq_lst = threed_data(uc, f1, f2)
        draw_qqplot(xq_lst, yq_lst)

    if input_map1.endswith((".cif")):
        fobs1, resol1 = pass_cif_sf(input_map1)
        fobs2, resol2 = pass_cif_sf(input_map2)
        xq_lst, yq_lst = make_xycols(fobs1, fobs2, resol1, resol2)
        draw_qqplot(xq_lst, yq_lst)
