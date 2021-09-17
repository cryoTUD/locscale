"""
Author: "Rangana Warshamanage, Garib N. Murshudov"
MRC Laboratory of Molecular Biology
    
This software is released under the
Mozilla Public License, version 2.0; see LICENSE.
"""
from __future__ import absolute_import, division, print_function, unicode_literals
import numpy as np

# Fullmap from halfmaps


def half2full(ar1, ar2):
    assert ar1.shape == ar2.shape
    hf1 = np.fft.fftn(ar1)
    hf2 = np.fft.fftn(ar2)
    return np.real(np.fft.ifftn((hf1 + hf2) / 2.0))

