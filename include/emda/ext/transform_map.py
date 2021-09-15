"""
Author: "Rangana Warshamanage, Garib N. Murshudov"
MRC Laboratory of Molecular Biology
    
This software is released under the
Mozilla Public License, version 2.0; see LICENSE.
"""

from __future__ import absolute_import, division, print_function, unicode_literals
import numpy as np
from emda import core, ext
import fcodes_fast


def write_mrc2(mapdata, filename, unit_cell, map_origin):
    import mrcfile as mrc

    data2write = np.real(np.fft.ifftshift(np.fft.ifftn(np.fft.ifftshift(mapdata))))
    # removing outer regions
    nx, ny, nz = data2write.shape
    assert nx == ny == nz
    dx = int(nx / 4)
    xx = int(nx / 2)
    newdata = data2write[dx : dx + xx, dx : dx + xx, dx : dx + xx]
    file = mrc.new(
        name=filename, data=np.float32(newdata), compression=None, overwrite=True
    )
    file.header.cella.x = unit_cell[0]
    file.header.cella.y = unit_cell[1]
    file.header.cella.z = unit_cell[2]
    file.header.nxstart = map_origin[0]
    file.header.nystart = map_origin[1]
    file.header.nzstart = map_origin[2]
    file.close()


def map_transform(mapname, t, r, ax, outname="transformed.mrc", mode="real"):
    from emda.ext.mapfit import utils
    from scipy.ndimage.interpolation import shift

    uc, arr, origin = core.iotools.read_map(mapname)
    if mode == "fourier":
        # in Fourier space - EMDA
        theta = [ax, r]
        q = core.quaternions.get_quaternion(theta)
        q = q / np.sqrt(np.dot(q, q))
        rotmat = core.quaternions.get_RM(q)
        rotmat = np.transpose(rotmat)
        print("rotation matrix:", rotmat)
        print(
            "Applied rotation in degrees [Euler angles]: ",
            core.quaternions.rotationMatrixToEulerAngles(rotmat) * 180.0 / np.pi,
        )
        print(
            "Applied rotation in degrees [Overall]: ",
            np.arccos((np.trace(rotmat) - 1) / 2) * 180.0 / np.pi,
        )
        print("Applied translation in Angstrom: ", t)
        hf = np.fft.fftshift(np.fft.fftn(arr))
        cx, cy, cz = hf.shape
        if np.sqrt(np.dot(t, t)) < 1.0e-3:
            st = 1.0
        else:
            t = (np.asarray(t)) / uc[:3]
            st, _, _, _ = fcodes_fast.get_st(cx, cy, cz, t)
        print("Applied translation in Angstrom: ", t * uc[:3])
        transformed_map = np.real(
            np.fft.ifftn(
                np.fft.ifftshift(utils.get_FRS(uc, rotmat, hf * st)[:, :, :, 0])
            )
        )

    if mode == "real":
        # real space map interpolation - fcodes_fast
        t_frac = utils.ang2frac(t, uc[0] / arr.shape[0])
        theta = [ax, r]
        q = core.quaternions.get_quaternion(theta)
        q = q / np.sqrt(np.dot(q, q))
        rotmat = core.quaternions.get_RM(q)
        print("rotation matrix:", rotmat)
        print(
            "Applied rotation in degrees [Euler angles]: ",
            core.quaternions.rotationMatrixToEulerAngles(rotmat) * 180.0 / np.pi,
        )
        print(
            "Applied rotation in degrees [Overall]: ",
            np.arccos((np.trace(rotmat) - 1) / 2) * 180.0 / np.pi,
        )
        print("Applied translation in Angstrom: ", t)
        nx, ny, nz = arr.shape
        transformed_map = fcodes_fast.trilinear_map(rotmat, arr, nx, ny, nz)
        transformed_map = shift(transformed_map, t_frac)

    """# using scipy - rotate function
    from scipy import ndimage
    f = np.fft.fftshift(np.fft.fftn(arr))
    frs_real = ndimage.rotate(f.real, r, axes=(1,2), reshape=False)
    frs_imag = ndimage.rotate(f.imag, r, axes=(1,2), reshape=False)
    frs = frs_real + 1j*frs_imag
    transformed_map = np.real(np.fft.ifftn(np.fft.ifftshift(frs)))
    angle_rad = np.deg2rad(r)
    c, s = np.cos(angle_rad), np.sin(angle_rad)
    rot_matrix = np.array([[c, s],
                              [-s, c]])
    print('rotation matrix:', rot_matrix)"""

    core.iotools.write_mrc(transformed_map, outname, uc, origin)
    return transformed_map


def map_transform_using_rotmat(arr, t, rotmat):
    # rotmat output by EMDA overlay and final translation(in fractions) are
    # directly used.
    from emda.ext.mapfit import utils

    fmap = np.fft.fftshift(np.fft.fftn(arr))
    nx, ny, nz = fmap.shape
    st, _, _, _ = fcodes_fast.get_st(nx, ny, nz, t)
    transformed_map = np.real(
        np.fft.ifftn(
            np.fft.ifftshift(
                utils.get_FRS(rotmat, fmap * st, interp="cubic")[:, :, :, 0]
            )
        )
    )
    return transformed_map
