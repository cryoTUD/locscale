"""
Author: "Rangana Warshamanage, Garib N. Murshudov"
MRC Laboratory of Molecular Biology
    
This software is released under the
Mozilla Public License, version 2.0; see LICENSE.
"""

from __future__ import absolute_import, division, print_function, unicode_literals
import numpy as np
import mrcfile as mrc


def read_mrc(mapnamelist, masknamelist):
    nmaps = len(mapnamelist)
    maplist = []
    masklist = []
    for i in range(nmaps):
        fid = mrc.open(mapnamelist[i])
        mid = mrc.open(masknamelist[i])
        if i == 0:
            map_orig = [
                1 * fid.header.nxstart,
                1 * fid.header.nystart,
                1 * fid.header.nzstart,
            ]
            cell = np.array(fid.header.cella)
            a = cell["x"]
            b = cell["y"]
            c = cell["z"]
            uc = np.asarray([a * 1.0, b * 1.0, c * 1.0, 90.0, 90.0, 90.0])
        mapar = np.asarray(fid.data, dtype="float")
        maskar = np.asarray(mid.data, dtype="float")
        mapar_centered = np.fft.fftshift(mapar)
        maskar_centered = np.fft.fftshift(maskar)
        maplist.append(mapar_centered)
        masklist.append(maskar_centered)
        fid.close()
        mid.close()
    return maplist, masklist, uc, map_orig


def write_mrc(filename, mapdata, unit_cell, map_origin):
    data2write = np.fft.ifftshift(mapdata)
    file = mrc.new(
        name=filename, data=np.float32(data2write), compression=None, overwrite=True
    )
    file.header.cella.x = unit_cell[0]
    file.header.cella.y = unit_cell[1]
    file.header.cella.z = unit_cell[2]
    file.header.nxstart = map_origin[0]
    file.header.nystart = map_origin[1]
    file.header.nzstart = map_origin[2]
    file.close()


def option_one(maplist, masklist, uc, map_orig):
    import numpy.ma as ma

    nx, ny, nz = maplist[0].shape
    w = np.zeros(shape=(nx, ny, nz), dtype="float")
    sum_MaskMap = np.zeros(shape=(nx, ny, nz), dtype="float")
    for imask in masklist:
        w += imask
    w_ma = ma.masked_less_equal(w, 0.0)
    nmaps = len(maplist)
    for i in range(nmaps):
        sum_MaskMap += masklist[i] * maplist[i]
    MapAver = sum_MaskMap / w_ma
    # Write out MapAver_filled into .mrc
    write_mrc("AverageMask_1.mrc", MapAver.filled(0.0), uc, map_orig)


def option_two(maplist, masklist, uc, map_orig):
    nx, ny, nz = maplist[0].shape
    maxMaskMap = np.zeros(shape=(nx, ny, nz), dtype="float")
    nmaps = len(maplist)
    for i in range(nmaps):
        if i == 0:
            maxMaskMap = masklist[i] * maplist[i]
        else:
            maxMaskMap = np.maximum(maxMaskMap, masklist[i] * maplist[i])
    # Write out MapAver_filled into .mrc
    write_mrc("MaxMaskMap_2.mrc", maxMaskMap, uc, map_orig)


def option_three(maplist, masklist, uc, map_orig):
    import numpy.ma as ma

    nx, ny, nz = maplist[0].shape
    sum_MaskMap = np.zeros(shape=(nx, ny, nz), dtype="float")
    sum_maxMask = np.zeros(shape=(nx, ny, nz), dtype="float")
    nmaps = len(maplist)
    for i in range(nmaps):
        sum_MaskMap += masklist[i] * maplist[i]
        sum_maxMask += masklist[i] * np.amax(masklist[i])
    sum_maxMask_ma = ma.masked_less_equal(sum_maxMask, 0.0)
    MapAver = sum_MaskMap / sum_maxMask_ma
    # Write out MapAver_filled into .mrc
    write_mrc("AverageMaxMask_3.mrc", MapAver.filled(0.0), uc, map_orig)


def main(mapslist, masklist):
    assert len(mapslist) == len(masklist)
    maps, masks, uc, map_orig = read_mrc(mapslist, masklist)
    print("Please wait until maps are ready...")
    option_one(maps, masks, uc, map_orig)
    option_two(maps, masks, uc, map_orig)
    option_three(maps, masks, uc, map_orig)
    print("Maps were written!!")
