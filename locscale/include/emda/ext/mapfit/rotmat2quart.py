"""
Author: "Rangana Warshamanage, Garib N. Murshudov"
MRC Laboratory of Molecular Biology
    
This software is released under the
Mozilla Public License, version 2.0; see LICENSE.
"""

# rotation matrix to quaternion
# http://www.euclideanspace.com/maths/geometry/rotations/conversions/matrixToQuaternion/index.htm

# assumes the upper 3x3 of m is a pure rotation matrix (i.e, unscaled)

from __future__ import absolute_import, division, print_function, unicode_literals
import math
import numpy as np
from emda import core


def rot2quart(rm):

    m11 = rm[0, 0]
    m12 = rm[0, 1]
    m13 = rm[0, 2]
    m21 = rm[1, 0]
    m22 = rm[1, 1]
    m23 = rm[1, 2]
    m31 = rm[2, 0]
    m32 = rm[2, 1]
    m33 = rm[2, 2]

    trace = m11 + m22 + m33

    q = np.zeros(4, dtype=np.float64)

    if trace > 0:
        s = 0.5 / math.sqrt(trace + 1.0)
        q[0] = 0.25 / s
        q[1] = (m32 - m23) * s
        q[2] = (m13 - m31) * s
        q[3] = (m21 - m12) * s
    elif m11 > m22 and m11 > m33:
        s = 2.0 * math.sqrt(1.0 + m11 - m22 - m33)
        q[0] = (m32 - m23) / s
        q[1] = 0.25 * s
        q[2] = (m12 + m21) / s
        q[3] = (m13 + m31) / s
    elif m22 > m33:
        s = 2.0 * math.sqrt(1.0 + m22 - m11 - m33)
        q[0] = (m13 - m31) / s
        q[1] = (m12 + m21) / s
        q[2] = 0.25 * s
        q[3] = (m23 + m32) / s
    else:
        s = 2.0 * math.sqrt(1.0 + m33 - m11 - m22)
        q[0] = (m21 - m12) / s
        q[1] = (m13 + m31) / s
        q[2] = (m23 + m32) / s
        q[3] = 0.25 * s
    return q


if __name__ == "__main__":
    # q_init = np.array([1.0, 0.0, 0.0, 0.0], dtype=np.float64)
    theta = [(0.5, 0.5, 0), 180.0]
    q_init = core.quaternions.get_quaternion(theta)
    print(q_init)
    q = rot2quart(core.quaternions.get_RM(q_init))
    print(q)
