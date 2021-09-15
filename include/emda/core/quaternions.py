import numpy as np
import math


def get_quaternion(theta):
    # angle should be in degrees
    rv = theta[0]
    rv = np.asarray(rv, dtype='float')
    rv = rv / np.sqrt(np.dot(rv, rv))
    ang = math.sin(np.deg2rad(theta[1]) / 2.0)
    q1, q2, q3 = rv[0] * ang, rv[1] * ang, rv[2] * ang
    q0 = math.cos(np.deg2rad(theta[1]) / 2.0)
    # Constraint to creat quaternion
    #q0 = np.sqrt(1 - q1 * q1 - q2 * q2 - q3 * q3)
    # Quaternion
    q = np.array([q0, q1, q2, q3], dtype=np.float64)
    return q


def quart2axis(q):
    q = q / np.sqrt(np.dot(q, q))
    angle = 2 * math.acos(q[0])
    s = math.sqrt(1.0 - q[0] * q[0])
    if s < 1e-3:
        x, y, z = q[1], q[2], q[3]
    else:
        x, y, z = q[1] / s, q[2] / s, q[3] / s
    return [x, y, z, angle]


def q_normalised(rv):
    q02 = 1 - np.sum(n * n for n in rv[1:])
    q2 = np.array([q02, abs(rv[1]), abs(rv[2]), abs(rv[3])], dtype=np.float64)
    q = np.sqrt(q2)
    return q


def get_RM(q):
    RM = np.array(
        [
            [
                1.0 - 2.0 * (q[2] ** 2 + q[3] ** 2),
                2.0 * (q[1] * q[2] - q[3] * q[0]),
                2.0 * (q[1] * q[3] + q[2] * q[0]),
            ],
            [
                2.0 * (q[1] * q[2] + q[3] * q[0]),
                1.0 - 2.0 * (q[1] ** 2 + q[3] ** 2),
                2.0 * (q[2] * q[3] - q[1] * q[0]),
            ],
            [
                2.0 * (q[1] * q[3] - q[2] * q[0]),
                2.0 * (q[2] * q[3] + q[1] * q[0]),
                1.0 - 2.0 * (q[1] ** 2 + q[2] ** 2),
            ],
        ],
        dtype=np.float64,
    )
    return RM


def derivatives_wrt_q(q):
    # First-order derivatives of R wrt q
    dRdq1 = np.array(
        [
            [
                0.0,
                2.0 * q[2] + 2.0 * (q[1] * q[3]) / q[0],
                2.0 * q[3] - 2.0 * (q[1] * q[2]) / q[0],
            ],
            [
                2.0 * q[2] - 2.0 * (q[1] * q[3]) / q[0],
                -2.0 * q[1],
                -2.0 * q[0] + 2.0 * (q[1] * q[1]) / q[0],
            ],
            [
                2.0 * q[3] + 2.0 * (q[1] * q[2]) / q[0],
                2.0 * q[0] - 2.0 * (q[1] * q[1]) / q[0],
                -2.0 * q[1],
            ],
        ],
        dtype=np.float64,
    )
    dRdq2 = np.array(
        [
            [
                -2.0 * q[2],
                2.0 * q[1] + 2.0 * (q[2] * q[3]) / q[0],
                2.0 * q[0] - 2.0 * (q[2] * q[2]) / q[0],
            ],
            [
                2.0 * q[1] - 2.0 * (q[2] * q[3]) / q[0],
                0.0,
                2.0 * q[3] - 2.0 * (q[1] * q[2]) / q[0],
            ],
            [
                -2.0 * q[0] + 2.0 * (q[2] * q[2]) / q[0],
                2.0 * q[3] - 2.0 * (q[1] * q[2]) / q[0],
                -2.0 * q[1],
            ],
        ],
        dtype=np.float64,
    )
    dRdq3 = np.array(
        [
            [
                -2.0 * q[3],
                -2.0 * q[0] + 2.0 * (q[3] * q[3]) / q[0],
                2.0 * q[1] - 2.0 * (q[2] * q[3]) / q[0],
            ],
            [
                2.0 * q[0] - 2.0 * (q[3] * q[3]) / q[0],
                -2.0 * q[3],
                2.0 * q[2] - 2.0 * (q[1] * q[3]) / q[0],
            ],
            [
                2.0 * q[1] + 2.0 * (q[2] * q[3]) / q[0],
                2.0 * q[2] - 2.0 * (q[1] * q[3]) / q[0],
                0.0,
            ],
        ],
        dtype=np.float64,
    )
    return np.array([dRdq1, dRdq2, dRdq3])


# source https://www.learnopencv.com/rotation-matrix-to-euler-angles/
# Checks if a matrix is a valid rotation matrix.
def isRotationMatrix(R):
    Rt = np.transpose(R)
    shouldBeIdentity = np.dot(Rt, R)
    I = np.identity(3, dtype=R.dtype)
    n = np.linalg.norm(I - shouldBeIdentity)
    return n < 1e-6


# Calculates rotation matrix to euler angles
# The result is the same as MATLAB except the order
# of the euler angles ( x and z are swapped ).
def rotationMatrixToEulerAngles(R):
    assert isRotationMatrix(R)
    sy = math.sqrt(R[0, 0] * R[0, 0] + R[1, 0] * R[1, 0])
    singular = sy < 1e-6
    if not singular:
        x = math.atan2(R[2, 1], R[2, 2])
        y = math.atan2(-R[2, 0], sy)
        z = math.atan2(R[1, 0], R[0, 0])
    else:
        x = math.atan2(-R[1, 2], R[1, 1])
        y = math.atan2(-R[2, 0], sy)
        z = 0
    return np.array([x, y, z])


def rotmat_from_axisangle(axis, theta):
    """
    Returns rotation matrix of counterclockwise rotation about axis.

    Rotation matrix is obtained from expanding Rodiguez formula.
    This requires angle in radians.
    """
    axis = np.asarray(axis)
    axis = axis / math.sqrt(np.dot(axis, axis))
    x, y, z = axis
    xx, yy, zz = x * x, y * y, z * z
    xy, xz, yz = x * y, x * z, y * z
    ca, sa = math.cos(theta), math.sin(theta)
    aa = 1 - ca
    return np.array(
        [
            [ca + xx * aa, xy * aa - z * sa, xz * aa + y * sa],
            [xy * aa + z * sa, ca + yy * aa, yz * aa - x * sa],
            [xz * aa - y * sa, yz * aa + x * sa, ca + zz * aa],
        ],
        dtype="float",
    )


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


def quaternion_inv(q):
    q_conjg = np.array([q[0], -q[1], -q[2], -q[3]], dtype='float')
    q_inv = q_conjg / np.dot(q, q)
    return q_inv


def quaternion_multiply(quaternion1, quaternion0):
    #https://stackoverflow.com/questions/39000758/how-to-multiply-two-quaternions-by-python-or-numpy
    w0, x0, y0, z0 = quaternion0
    w1, x1, y1, z1 = quaternion1
    return np.array([-x1 * x0 - y1 * y0 - z1 * z0 + w1 * w0,
                     x1 * w0 + y1 * z0 - z1 * y0 + w1 * x0,
                     -x1 * z0 + y1 * w0 + z1 * x0 + w1 * y0,
                     x1 * y0 - y1 * x0 + z1 * w0 + w1 * z0], dtype=np.float64)