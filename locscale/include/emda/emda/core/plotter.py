"""
Author: "Rangana Warshamanage, Garib N. Murshudov"
MRC Laboratory of Molecular Biology

This software is released under the
Mozilla Public License, version 2.0; see LICENSE.
"""

from __future__ import absolute_import, division, print_function, unicode_literals
import numpy as np
import matplotlib

#matplotlib.use("Agg")
matplotlib.use(matplotlib.get_backend())
from matplotlib import pyplot as plt


def plot_line_log(res_arr, arr):
    d2vec = []
    for i in range(len(res_arr)):
        d2vec.append(float(1 / (res_arr[i])))
    fig, ax = plt.subplots(1, 1)
    plt.plot(d2vec, np.log10(arr), "r")
    xmin = np.min(d2vec)
    xmax = np.max(d2vec)
    ax.set_xlim(xmin, xmax)
    ax.set_xticks(d2vec)
    str1 = []
    for e in res_arr:
        str1.append("{:.3}".format(e))
    for n in range(len(str1)):
        if n % 5 != 0:
            str1[n] = " "
    ax.set_xticklabels(str1, rotation=90, fontsize=9)
    plt.xlabel("Resolution")
    plt.ylabel("log.Variance(Fo)")
    plt.show()


def plot_nlines_log(
    res_arr,
    list_arr,
    # curve_label=["Signal","Noise","Predicted"],
    curve_label=None,
    mapname=None,
    plot_title=None,
    xlabel="Resolution ($\AA$)",
    ylabel="log.Variance(Fo)",
):
    if mapname is None:
        mapname = "variances.eps"
    bin_arr = np.arange(len(res_arr))
    fig = plt.figure(figsize=(6, 4))
    ax1 = fig.add_subplot(111)
    for icurve in range(len(list_arr)):
        data = list_arr[icurve]
        if curve_label is None:
            label = "Set#" + str(icurve)
            ax1.plot(bin_arr, np.log10(data), label=label, linewidth=2)
        else:
            ax1.plot(bin_arr, np.log10(data), label=curve_label[icurve], linewidth=2)
    pos = np.array(ax1.get_xticks(), dtype=np.int)
    n_bins = res_arr.shape[0]
    pos[pos < 0] = 0
    pos[pos >= n_bins] = n_bins - 1
    ax1.set_xticklabels(np.round(res_arr[pos], decimals=2))
    ax1.set_xlabel(xlabel)
    plt.ylabel(ylabel)
    #ax1.set_ylim([-2, 6])
    font = {"family": "serif", "color": "black", "weight": "bold", "size": 16}
    plt.legend(loc=0)
    font = {"family": "serif", "color": "black", "weight": "bold", "size": 12}
    if plot_title is not None:
        plt.title(plot_title, fontdict=font)
    plt.savefig(mapname, format="eps", dpi=300)
    #plt.show()
    plt.close()


def plot_line(mapname, res_arr, arr):
    bin_arr = np.arange(len(res_arr))
    fig = plt.figure(figsize=(6, 4))
    ax1 = fig.add_subplot(111)
    ax1.plot(bin_arr, arr, linewidth=2)
    xmin = np.min(bin_arr)
    xmax = np.max(bin_arr)
    # plt.plot((xmin, xmax), (0.143, 0.143), 'k--')
    pos = np.array(ax1.get_xticks(), dtype=np.int)
    n_bins = res_arr.shape[0]
    pos[pos < 0] = 0
    pos[pos >= n_bins] = n_bins - 1
    ax1.set_xticklabels(np.round(res_arr[pos], decimals=2))
    ax1.set_xlabel("Resolution ($\AA$)")
    plt.ylabel("Fourier Shell Correlation")
    font = {"family": "serif", "color": "black", "weight": "bold", "size": 16}
    plt.savefig(mapname, format="eps", dpi=300)
    plt.show()


def plot_nlines(
    res_arr,
    list_arr,
    mapname="halfmap_fsc.eps",
    curve_label=None,
    fscline=0.143,
    plot_title=None,
):
    # from mpl_toolkits.axes_grid1 import host_subplot
    # import mpl_toolkits.axisartist as AA

    if curve_label is None:
        curve_label = ["halfmap_fsc"]
    bin_arr = np.arange(len(res_arr))
    fig = plt.figure(figsize=(6, 4))
    # ax1 = host_subplot(111,axes_class=AA.Axes)
    ax1 = fig.add_subplot(111)
    for icurve in range(len(list_arr)):
        ax1.plot(bin_arr, list_arr[icurve], label=curve_label[icurve], linewidth=2)
    xmin = np.min(bin_arr)
    xmax = np.max(bin_arr)
    plt.plot(
        (xmin, xmax), (float(fscline), float(fscline)), color="gray", linestyle=":"
    )
    plt.plot((xmin, xmax), (0.0, 0.0), color="black", linestyle=":")
    pos = np.array(ax1.get_xticks(), dtype=np.int)
    n_bins = res_arr.shape[0]
    pos[pos < 0] = 0
    pos[pos >= n_bins] = n_bins - 1
    #
    import matplotlib.ticker as mticker
    label_format = '{:5.2f}'
    ax1.xaxis.set_major_locator(mticker.FixedLocator(pos))
    ax1.set_xticklabels([label_format.format(res_arr[x]) for x in pos])
    #
    #ax1.set_xticklabels(np.round(res_arr[pos], decimals=2))
    ax1.set_xlabel("Resolution ($\AA$)")
    plt.legend(loc=0)
    plt.ylabel("Fourier Shell Correlation")
    font = {"family": "serif", "color": "black", "weight": "bold", "size": 12}
    if plot_title is not None:
        plt.title(plot_title, fontdict=font)
    plt.savefig(mapname, format="eps", dpi=300)
    plt.close()


def plot_nlines2(
    res_arr, list_arr, mapname="halfmap_fsc.eps", curve_label=None
):
    # from mpl_toolkits.axes_grid1 import host_subplot
    # import mpl_toolkits.axisartist as AA

    if curve_label is None:
        curve_label = ["halfmap_fsc"]
    bin_arr = np.arange(len(res_arr))
    fig = plt.figure(figsize=(6, 4))
    # ax1 = host_subplot(111,axes_class=AA.Axes)
    ax1 = fig.add_subplot(111)
    for icurve in range(len(list_arr)):
        ax1.plot(bin_arr, list_arr[icurve], label=curve_label[icurve], linewidth=2)
    xmin = np.min(bin_arr)
    xmax = np.max(bin_arr)
    plt.plot((xmin, xmax), (0.5, 0.5), "k--")

    pos = np.array(ax1.get_xticks(), dtype=np.int)
    n_bins = res_arr.shape[0]
    pos[pos < 0] = 0
    pos[pos >= n_bins] = n_bins - 1
    ax1.set_xticklabels(np.round(res_arr[pos], decimals=2))
    ax1.set_xlabel("Resolution (1/$\AA$)")
    plt.legend(loc=0)
    plt.ylabel("Fourier Shell Correlation")
    font = {"family": "serif", "color": "black", "weight": "bold", "size": 16}
    plt.savefig(mapname, format="eps", dpi=300)


def plot_3d(data):
    # from mpl_toolkits.mplot3d import Axes3D  # implicit import

    fig = plt.figure()
    ax = fig.gca(projection="3d")
    # nan_mask = Z > 0
    r = len(data)
    x_dim = r
    y_dim = r
    X = np.linspace(1, x_dim, x_dim)
    Y = np.linspace(1, y_dim, y_dim)
    ax.scatter(X, Y, data, marker=".", color="r", s=3)
    plt.xlabel("X")
    plt.ylabel("Y")
    ax.set_zlabel("Z")
    ax.axis("equal")
    ax.axis("tight")
    plt.show()


def plot3dnonzero(data3d):
    # from mpl_toolkits.mplot3d import Axes3D

    dim = data3d.shape[0]
    z, x, y = data3d.nonzero()
    fig = plt.figure()
    ax = fig.add_subplot(111, projection="3d")
    ax.scatter(x, y, -z, zdir="z", c="red", marker=".")
    plt.show()


def plot4dcolor(data):
    # from mpl_toolkits.mplot3d import Axes3D

    ih, ik, il = data.shape
    h = np.arange(-ih / 2, ih / 2)
    k = np.arange(-ik / 2, ik / 2)
    l = np.arange(-il / 2, il / 2)
    y, x, z = np.meshgrid(k, h, l)
    c = data.flatten()
    fig = plt.figure()
    ax = fig.add_subplot(111, projection="3d")
    ax.scatter(x, y, z, c=c, cmap=plt.hot())
    # ax.scatter(x, y, z, c, color='b')
    plt.show()


def contour_nplot(nx, ny, lst_data):
    plt.style.use("seaborn-white")
    plt.figure(figsize=(30, 10))
    x = np.linspace(0, nx, nx)
    y = np.linspace(0, ny, ny)
    X, Y = np.meshgrid(x, y)
    lst_minval = []
    lst_maxval = []
    for imap in lst_data:
        lst_minval.append(np.amin(imap))
        lst_maxval.append(np.amax(imap))
    min_val = min(lst_minval)
    max_val = max(lst_maxval)
    i = 0
    n = len(lst_data)
    # fig = plt.figure(figsize=(30,10))
    from mpl_toolkits.axes_grid1 import make_axes_locatable

    fig, (ax1, ax2, ax3) = plt.subplots(ncols=3)
    for imap in lst_data:
        i = i + 1
        plt.subplot(1, n, i)
        if i == 3:
            img1 = plt.contourf(
                X, Y, imap, 200, cmap="RdGy", levels=np.linspace(min_val, max_val, 20)
            )
            divider = make_axes_locatable(ax3)
            cax1 = divider.append_axes("right", size="5%", pad=0.05)
            fig.colorbar(img1, cax=cax1)
        else:
            plt.contourf(
                X, Y, imap, 200, cmap="RdGy", levels=np.linspace(min_val, max_val, 20)
            )
    plt.show()


def contour_nplot2(data_3d):
    plt.style.use("seaborn-white")
    # plt.figure(figsize=(10,10))
    nx, ny, nz = data_3d.shape

    x = np.linspace(0, nx, nx)
    y = np.linspace(0, ny, ny)
    X, Y = np.meshgrid(x, y)
    min_val = np.amin(data_3d)
    max_val = np.amax(data_3d)
    # from mpl_toolkits.axes_grid1 import make_axes_locatable
    # fig, (ax1, ax2, ax3) = plt.subplots(ncols=1)
    for i in range(data_3d.shape[2]):
        plt.title("Z= %i" % i)
        plt.contourf(
            X,
            Y,
            data_3d[:, :, i],
            200,
            cmap="RdGy",
            levels=np.linspace(min_val, max_val, 20),
        )
        plt.axis("equal")
        plt.show()


def plot_from_twofiles_csv(filename, labels):
    import pandas

    data_refmac = pandas.read_csv(filename[0], delimiter=" ")
    data_other = pandas.read_csv(filename[1], delimiter=" ")
    bf_arr = abs(data_other[labels[0]])
    blur_cc_other = data_other[labels[1]]
    blur_cc_refmac = data_refmac[labels[1]]
    # bf_arr = np.array([0,10,20,30,40,50,60,70,80,90,100,120,140,160])
    fig = plt.figure(figsize=(6, 4))
    ax1 = fig.add_subplot(111)
    ax1.plot(bf_arr, blur_cc_other, label="Other", linewidth=2)
    ax1.plot(bf_arr, blur_cc_refmac, label="Refmac", linewidth=2)
    xmin = np.min(bf_arr)
    xmax = np.max(bf_arr)
    plt.plot((xmin, xmax), (1.0, 1.0), "k--")
    ax1.set_xlabel("B factor ($\AA^2$)")
    plt.legend(loc=0)
    plt.ylabel("Overall Correlation Coefficient")
    plt.savefig("overall_CC_both.eps", format="eps", dpi=300)
    plt.show()


def plot_fval(fval_list,i=0):
    fig = plt.figure(figsize=(6, 4))
    ax1 = fig.add_subplot(111)
    ax1.plot(fval_list, label="fval", linewidth=2)
    plt.legend(loc=0)
    plt.ylabel("functional value")
    outname = "fval_mcycle_%i" %i
    plt.savefig(outname+".eps", format="eps", dpi=300)
    plt.show()    