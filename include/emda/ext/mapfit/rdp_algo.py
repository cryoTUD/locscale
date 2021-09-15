"""
Author: "Rangana Warshamanage, Garib N. Murshudov"
MRC Laboratory of Molecular Biology
    
This software is released under the
Mozilla Public License, version 2.0; see LICENSE.
"""


# The Ramer-Douglas-Peucker algorithm roughly ported from the pseudo-code provided
# by http://en.wikipedia.org/wiki/Ramer-Douglas-Peucker_algorithm
# Implementation: https://github.com/sebleier/RDP/blob/master/__init__.py
# MIT License

from __future__ import absolute_import, division, print_function, unicode_literals
from math import sqrt
import numpy as np


def run_rdp(x1d_np, y1d_np, epsilon):
    merged_list = tuple(zip(x1d_np.tolist(), y1d_np.tolist()))
    results = rdp(merged_list, epsilon)
    x_new, y_new = zip(*results)
    return np.asarray(x_new), np.asarray(y_new)


def distance(a, b):
    return sqrt((a[0] - b[0]) ** 2 + (a[1] - b[1]) ** 2)


def point_line_distance(point, start, end):
    if start == end:
        return distance(point, start)
    else:
        n = abs(
            (end[0] - start[0]) * (start[1] - point[1])
            - (start[0] - point[0]) * (end[1] - start[1])
        )
        d = sqrt((end[0] - start[0]) ** 2 + (end[1] - start[1]) ** 2)
        return n / d


def rdp(points, epsilon):
    """Reduces a series of points to a simplified version that loses detail, but
    maintains the general shape of the series.
    """
    dmax = 0.0
    index = 0
    for i in range(1, len(points) - 1):
        d = point_line_distance(points[i], points[0], points[-1])
        if d > dmax:
            index = i
            dmax = d

    if dmax >= epsilon:
        results = rdp(points[: index + 1], epsilon)[:-1] + rdp(points[index:], epsilon)
    else:
        results = [points[0], points[-1]]

    return results
