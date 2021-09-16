"""
Author: "Rangana Warshamanage, Garib N. Murshudov"
MRC Laboratory of Molecular Biology
    
This software is released under the
Mozilla Public License, version 2.0; see LICENSE.
"""

from __future__ import absolute_import, division, print_function, unicode_literals
import numpy as np
import math


def cc_overall_realsp(map1, map2, mask=None):
    # Overall correlation - real space
    if mask is not None:
        map1 = map1 * mask
        map2 = map2 * mask
        mean1 = np.sum(map1) / np.sum(mask)
        mean2 = np.sum(map2) / np.sum(mask)
        mean1 = mean1 * mask
        mean2 = mean2 * mask
    else:
        mean1 = np.mean(map1)
        mean2 = np.mean(map2)
    covar = np.sum(map1 * map2 - (mean1 * mean2))
    var1 = np.sum(map1 * map1 - (mean1 * mean1))
    var2 = np.sum(map2 * map2 - (mean2 * mean2))
    occ = covar / math.sqrt(var1 * var2)
    return 2.0 * occ / (1.0 + occ), occ


def cc_overall_fouriersp(f1, f2):
    # Overall correlation - fourier space
    f1_mean = np.mean(f1)
    f2_mean = np.mean(f2)
    covar = np.sum(f1 * np.conjugate(f2) - f1_mean * np.conjugate(f2_mean))
    var1 = np.sum(f1 * np.conjugate(f1) - f1_mean * np.conjugate(f1_mean))
    var2 = np.sum(f2 * np.conjugate(f2) - f2_mean * np.conjugate(f2_mean))
    occ = covar.real / math.sqrt(var1.real * var2.real)
    return 2.0 * occ / (1.0 + occ), occ
