"""
Author: "Rangana Warshamanage, Garib N. Murshudov"
MRC Laboratory of Molecular Biology

This software is released under the
Mozilla Public License, version 2.0; see LICENSE.
"""

from __future__ import absolute_import, division, print_function, unicode_literals
import os
import numpy as np
import fcodes_fast
import emda.ext.mapfit.utils as utils
from emda import core, ext
import emda.emda_methods as em
from numpy.fft import fftn, fftshift


# Overlay of several maps using maps. not use halfmaps
""" class EmmapOverlay:
    def __init__(self, hfmap_list, mask_list=None):
        self.hfmap_list = hfmap_list
        self.mask_list = mask_list
        self.map_unit_cell = None
        self.map_origin = None
        self.map_dim = None
        self.arr_lst = []
        self.ceo_lst = None
        self.cfo_lst = None
        self.cbin_idx = None
        self.cdim = None
        self.cbin = None
        self.com = False
        self.com1 = None
        self.box_centr = None
        self.fhf_lst = None
        self.nbin = None
        self.res_arr = None
        self.bin_idx = None
        self.fo_lst = None
        self.eo_lst = None
        self.totalvar_lst = None

    def load_maps(self, fobj):
        from scipy import ndimage
        from scipy.ndimage.interpolation import shift

        com = self.com
        cmask = False
        fhf_lst = []
        if self.mask_list is not None:
            if len(self.hfmap_list) != len(self.mask_list):
                print("map_list and mask_list must have the same size!")
                print("exiting program...")
                exit()
            for i in range(len(self.mask_list)):
                _, mask, _ = core.iotools.read_map(self.mask_list[i])
                mask = utils.set_dim_even(mask)
                uc, arr, origin = core.iotools.read_map(self.hfmap_list[i])
                arr = utils.set_dim_even(arr)
                try:
                    assert arr.shape == mask.shape
                except AssertionError:
                    print("Map and Mask Dimension mismatched!")
                    exit()
                arr = arr * mask
                if i == 0:
                    nx, ny, nz = arr.shape
                    map_origin = origin
                    uc_target = uc
                    target_dim = arr.shape
                    target_pix_size = uc_target[0] / target_dim[0]
                    if cmask:
                        corner_mask = utils.remove_unwanted_corners(uc, target_dim)
                    else:
                        corner_mask = 1.0
                    if com:
                        # temp = arr * mask
                        com1 = ndimage.measurements.center_of_mass(arr * (arr >= 0.0))
                        print("COM: ", com1)
                        box_centr = (nx // 2, ny // 2, nz // 2)
                        print(box_centr)
                        self.com1 = com1
                        self.box_centr = box_centr
                        arr_mvd = shift(arr, np.subtract(box_centr, com1))
                        self.arr_lst.append(arr_mvd * corner_mask)
                        fhf_lst.append(
                            np.fft.fftshift(
                                np.fft.fftn(np.fft.fftshift(arr_mvd * corner_mask))
                            )
                        )
                    else:
                        fhf_lst.append(
                            np.fft.fftshift(
                                np.fft.fftn(np.fft.fftshift(arr * corner_mask))
                            )
                        )
                else:
                    curnt_pix_size = uc[0] / arr.shape[0]
                    arr = core.iotools.resample2staticmap(
                        curnt_pix=curnt_pix_size,
                        targt_pix=target_pix_size,
                        targt_dim=target_dim,
                        arr=arr,
                    )
                    if com:
                        com1 = ndimage.measurements.center_of_mass(arr * (arr >= 0.0))
                        print("COM: ", com1)
                        arr = shift(arr, np.subtract(box_centr, com1))
                    self.arr_lst.append(arr * corner_mask)
                    fhf_lst.append(
                        np.fft.fftshift(np.fft.fftn(np.fft.fftshift(arr * corner_mask)))
                    )
            self.map_origin = map_origin
            self.map_unit_cell = uc_target
            self.map_dim = target_dim
            self.fhf_lst = fhf_lst
        if self.mask_list is None:
            for i in range(len(self.hfmap_list)):
                uc, arr, origin = core.iotools.read_map(self.hfmap_list[i])
                arr = utils.set_dim_even(arr)
                print("origin: ", origin)
                if i == 0:
                    nx, ny, nz = arr.shape
                    map_origin = origin
                    uc_target = uc
                    target_dim = arr.shape
                    target_pix_size = uc_target[0] / target_dim[0]
                    if com:
                        com1 = ndimage.measurements.center_of_mass(arr * (arr >= 0.0))
                        print("COM before centering: ", com1)
                        box_centr = (nx // 2, ny // 2, nz // 2)
                        print("BOX center: ", box_centr)
                        self.com1 = com1
                        self.box_centr = box_centr
                        arr = shift(arr, np.subtract(box_centr, com1))
                        self.arr_lst.append(arr)
                        core.iotools.write_mrc(
                            arr, "static_centered.mrc", uc_target, map_origin
                        )
                        print(
                            "COM after centering: ",
                            ndimage.measurements.center_of_mass(arr * (arr >= 0.0)),
                        )
                    self.arr_lst.append(arr)
                    fhf_lst.append(np.fft.fftshift(np.fft.fftn(np.fft.fftshift(arr))))
                else:
                    curnt_pix_size = uc[0] / arr.shape[0]
                    arr = core.iotools.resample2staticmap(
                        curnt_pix=curnt_pix_size,
                        targt_pix=target_pix_size,
                        targt_dim=target_dim,
                        arr=arr,
                    )
                    if com:
                        com1 = ndimage.measurements.center_of_mass(arr * (arr >= 0.0))
                        print("COM: ", com1)
                        arr = shift(arr, np.subtract(box_centr, com1))
                        core.iotools.write_mrc(
                            arr, "moving_centered.mrc", uc_target, map_origin
                        )
                        print(
                            "COM after centering: ",
                            ndimage.measurements.center_of_mass(arr * (arr >= 0.0)),
                        )
                    self.arr_lst.append(arr)
                    fhf_lst.append(np.fft.fftshift(np.fft.fftn(np.fft.fftshift(arr))))
            self.map_origin = map_origin
            self.map_unit_cell = uc_target
            self.map_dim = target_dim
            self.fhf_lst = fhf_lst

    def load_models(self):
        import numpy as np
        import emda.emda_methods as em

        dim_list = []
        model_list = self.hfmap_list
        shifted_model_list = []
        for i, model in enumerate(model_list):
            shifted_model_list.append("shifted" + str(i) + ".cif")
            dim_list.append(
                em.get_dim(model=model, shiftmodel=shifted_model_list[i])
            )
        dim = max(dim_list)

        fhf_lst = []
        unit_cell = np.array([dim, dim, dim, 90.0, 90.0, 90.0], dtype='float')
        for i, model in enumerate(shifted_model_list):
            modelmap = em.model2map(
                model,
                dim=[dim, dim, dim],
                resol=4,
                cell=[float(dim), float(dim), float(dim), 90.0, 90.0, 90.0],
            )
            em.write_mrc(
                np.fft.fftshift(modelmap),
                # modelmap,
                "modelmap_reboxed%s.mrc" % (i),
                unit_cell=unit_cell,
            )
            fhf_lst.append(np.fft.fftshift(np.fft.fftn(modelmap)))

        self.map_origin = [0, 0, 0]
        self.map_unit_cell = unit_cell
        self.map_dim = [dim, dim, dim]
        self.fhf_lst = fhf_lst """

class EmmapOverlay:
    def __init__(self, map_list, modelres=5.0, com=False, mask_list=None):
        self.map_list = map_list
        self.mask_list = mask_list
        self.modelres = modelres
        self.map_unit_cell = None
        self.map_origin = None
        self.map_dim = None
        self.pixsize = None
        self.arr_lst = []
        self.ceo_lst = None
        self.cfo_lst = None
        self.cbin_idx = None
        self.cdim = None
        self.cbin = None
        self.com = com
        self.com1 = None
        self.comlist = []
        self.box_centr = None
        self.fhf_lst = None
        self.nbin = None
        self.res_arr = None
        self.bin_idx = None
        self.fo_lst = None
        self.eo_lst = None
        self.totalvar_lst = None

    def load_maps(self, fobj):
        from scipy import ndimage
        from scipy.ndimage.interpolation import shift

        com = self.com
        cmask = False
        fhf_lst = []
        if self.mask_list is not None:
            if len(self.map_list) != len(self.mask_list):
                raise SystemExit("map_list and mask_list must have the same size!")
            for i in range(len(self.mask_list)):
                if i == 0:
                    _, mask, _ = em.get_data(self.mask_list[i])
                    mask = utils.set_dim_even(mask)
                    uc, arr, origin = em.get_data(self.map_list[i])
                    arr = utils.set_dim_even(arr)
                    try:
                        assert arr.shape == mask.shape
                    except AssertionError:
                        raise SystemExit("Map and Mask Dimension mismatched!")
                    arr = arr * mask
                    nx, ny, nz = arr.shape
                    map_origin = origin
                    uc_target = uc
                    target_dim = arr.shape
                    target_pix_size = []
                    for i in range(3):
                        target_pix_size.append(round(uc_target[i] / target_dim[i], 5))
                    #target_pix_size = uc_target[0] / target_dim[0]
                    if cmask:
                        corner_mask = utils.remove_unwanted_corners(uc, target_dim)
                    else:
                        corner_mask = 1.0
                    if com:
                        com1 = ndimage.measurements.center_of_mass(arr * (arr >= 0.0))
                        print("Before centering COM: ", com1)
                        box_centr = (nx // 2, ny // 2, nz // 2)
                        self.com1, self.box_centr = com1, box_centr
                        self.comlist.append(com1)
                        arr_mvd = shift(arr, np.subtract(box_centr, com1))
                        com1 = ndimage.measurements.center_of_mass(arr_mvd * (arr_mvd >= 0.0))
                        print("After centering COM: ", com1)
                        self.arr_lst.append(arr_mvd * corner_mask)
                        fhf_lst.append(fftshift(fftn(fftshift(arr_mvd * corner_mask))))
                    else:
                        fhf_lst.append(fftshift(fftn(fftshift(arr * corner_mask))))
                else:
                    uc, arr, origin = em.get_data(
                        self.map_list[i],
                        resol=self.modelres,
                        dim=target_dim,
                        uc=uc_target,
                        maporigin=map_origin,
                    )
                    arr = utils.set_dim_even(arr)
                    print("origin: ", origin)
                    _, mask, _ = em.get_data(self.mask_list[i])
                    mask = utils.set_dim_even(mask)
                    try:
                        assert arr.shape == mask.shape
                    except AssertionError:
                        raise SystemExit("Map and Mask Dimension mismatched!")
                    arr = arr * mask
                    curnt_pix_size = []
                    for i in range(3):
                        curnt_pix_size.append(round(uc[i] / arr.shape[i], 5))
                    #curnt_pix_size = uc[0] / arr.shape[0]
                    arr = core.iotools.resample2staticmap(
                        curnt_pix=curnt_pix_size,
                        targt_pix=target_pix_size,
                        targt_dim=target_dim,
                        arr=arr,
                    )
                    if com:
                        com1 = ndimage.measurements.center_of_mass(arr * (arr >= 0.0))
                        print("Before centering COM: ", com1)
                        self.comlist.append(com1)
                        arr = shift(arr, np.subtract(box_centr, com1))
                        com1 = ndimage.measurements.center_of_mass(arr * (arr >= 0.0))
                        print("After centering COM: ", com1)
                    self.arr_lst.append(arr * corner_mask)
                    fhf_lst.append(fftshift(fftn(fftshift(arr * corner_mask))))
            self.pixsize = target_pix_size
            self.map_origin = map_origin
            self.map_unit_cell = uc_target
            self.map_dim = target_dim
            self.fhf_lst = fhf_lst
        if self.mask_list is None:
            """ uc, arr1, origin = em.get_data(self.map_list[0])
            uc, arr2, origin = em.get_data(self.map_list[1])
            f1 = (fftshift(fftn(fftshift(arr1))))
            f2 = (fftshift(fftn(fftshift(arr2))))
            #f1 = (fftshift(fftn((arr1))))
            #f2 = (fftshift(fftn((arr2))))
            nbin, res_arr, bin_idx = core.restools.get_resolution_array(uc, f1)
            # test FSC calculation between maps
            f1f2_fsc = core.fsc.anytwomaps_fsc_covariance(
                    f1, 
                    f2, 
                    bin_idx, 
                    nbin)[0]
            for i in range(nbin):
                print(i, res_arr[i], f1f2_fsc[i])
            exit() """            
            for i in range(len(self.map_list)):
                if i == 0:
                    uc, arr, origin = em.get_data(self.map_list[i])
                    arr = utils.set_dim_even(arr)
                    print("origin: ", origin)
                    nx, ny, nz = arr.shape
                    map_origin = origin
                    uc_target = uc
                    target_dim = arr.shape
                    target_pix_size = []
                    for i in range(3):
                        target_pix_size.append(round(uc_target[i] / target_dim[i], 5))
                    #target_pix_size = uc_target[0] / target_dim[0]
                    if com:
                        com1 = ndimage.measurements.center_of_mass(arr * (arr >= 0.0))
                        print("COM before centering: ", com1)
                        box_centr = (nx // 2, ny // 2, nz // 2)
                        print("BOX center: ", box_centr)
                        self.com1 = com1
                        self.box_centr = box_centr
                        self.comlist.append(com1)
                        arr = shift(arr, np.subtract(box_centr, com1))
                        self.arr_lst.append(arr)
                        core.iotools.write_mrc(
                            arr, "static_centered.mrc", uc_target, map_origin
                        )
                        print(
                            "COM after centering: ",
                            ndimage.measurements.center_of_mass(arr * (arr >= 0.0)),
                        )
                    self.arr_lst.append(arr)
                    fhf_lst.append(fftshift(fftn(fftshift(arr))))
                else:
                    uc, arr, origin = em.get_data(
                        self.map_list[i],
                        resol=self.modelres,
                        dim=target_dim,
                        uc=uc_target,
                        maporigin=map_origin,
                    )
                    arr = utils.set_dim_even(arr)
                    em.write_mrc(arr, 'modelmap'+str(i)+'.mrc', uc, origin)
                    print("origin: ", origin)
                    curnt_pix_size = []
                    for i in range(3):
                        curnt_pix_size.append(round(uc[i] / arr.shape[i], 5))
                    #curnt_pix_size = uc[0] / arr.shape[0]
                    arr = core.iotools.resample2staticmap(
                        curnt_pix=curnt_pix_size,
                        targt_pix=target_pix_size,
                        targt_dim=target_dim,
                        arr=arr,
                    )
                    print("target_dim", target_dim)
                    print("arr.shape: ", arr.shape)
                    assert np.all(target_dim == arr.shape)
                    if com:
                        com1 = ndimage.measurements.center_of_mass(arr * (arr >= 0.0))
                        print("COM: ", com1)
                        self.comlist.append(com1)
                        arr = shift(arr, np.subtract(box_centr, com1))
                        core.iotools.write_mrc(
                            arr, "moving_centered.mrc", uc_target, map_origin
                        )
                        print(
                            "COM after centering: ",
                            ndimage.measurements.center_of_mass(arr * (arr >= 0.0)),
                        )
                    self.arr_lst.append(arr)
                    fhf_lst.append(fftshift(fftn(fftshift(arr))))
            self.pixsize = target_pix_size
            self.map_origin = map_origin
            self.map_unit_cell = uc_target
            self.map_dim = target_dim
            self.fhf_lst = fhf_lst

    def load_models(self):
        self.com = False
        dim_list = []
        model_list = self.map_list
        shifted_model_list = []
        for i, model in enumerate(model_list):
            shifted_model_list.append("shifted" + str(i) + ".cif")
            dim_list.append(
                em.get_dim(model=model, shiftmodel=shifted_model_list[i])
            )
        dim = max(dim_list)

        fhf_lst = []
        unit_cell = np.array([dim, dim, dim, 90.0, 90.0, 90.0], dtype='float')
        for i, model in enumerate(shifted_model_list):
            modelmap = em.model2map(
                modelxyz=model,
                dim=[dim, dim, dim],
                resol=4,
                cell=[float(dim), float(dim), float(dim), 90.0, 90.0, 90.0],
            )
            em.write_mrc(
                np.fft.fftshift(modelmap),
                # modelmap,
                "modelmap_reboxed%s.mrc" % (i),
                unit_cell=unit_cell,
            )
            fhf_lst.append(fftshift(fftn(modelmap)))
        self.pixsize = unit_cell[0] / dim
        self.map_origin = [0, 0, 0]
        self.map_unit_cell = unit_cell
        self.map_dim = [dim, dim, dim]
        self.fhf_lst = fhf_lst


    def calc_fsc_from_maps(self, fobj):
        # function for only two maps fitting
        nmaps = len(self.fhf_lst)
        fFo_lst = []
        fEo_lst = []
        fBTV_lst = []

        self.nbin, self.res_arr, self.bin_idx = core.restools.get_resolution_array(
            self.map_unit_cell, self.fhf_lst[0]
        )
        """ # test FSC calculation between maps - problematic
        f1f2_fsc = core.fsc.anytwomaps_fsc_covariance(
                self.fhf_lst[0], 
                self.fhf_lst[1], 
                self.bin_idx, 
                self.nbin)[0]
        for i in range(self.nbin):
            print(i, self.res_arr[i], f1f2_fsc[i])
        exit() """
        #
        for i in range(nmaps):
            _, _, _, totalvar, fo, eo = core.fsc.halfmaps_fsc_variance(
                self.fhf_lst[i], self.fhf_lst[i], self.bin_idx, self.nbin
            )
            fFo_lst.append(fo)
            fEo_lst.append(eo)
            fBTV_lst.append(totalvar)
        #
        self.fo_lst = fFo_lst
        self.eo_lst = fEo_lst
        self.totalvar_lst = fBTV_lst


# Overlay of several maps using halfmaps
class Overlay:
    def __init__(self, hfmap_list, mask_list=None):
        self.hfmap_list = hfmap_list
        self.mask_list = mask_list
        self.map_unit_cell = None
        self.map_origin = None
        self.map_dim = None
        self.arr_lst = []
        self.ceo_lst = None
        self.cfo_lst = None
        self.cbin_idx = None
        self.cdim = None
        self.cbin = None
        self.com1 = None
        self.box_centr = None
        self.fhf_lst = None
        self.nbin = None
        self.res_arr = None
        self.bin_idx = None
        self.fo_lst = None
        self.eo_lst = None
        self.totalvar_lst = None

    def load_maps(self, fobj):
        from scipy import ndimage
        from scipy.ndimage.interpolation import shift

        fhf_lst = []
        com = False
        if self.mask_list is not None:
            if len(self.hfmap_list) // 2 != len(self.mask_list):
                print(
                    "mask_list size is not equal to half the size of map_list!"
                )
                raise SystemExit("exiting program...")
            for i in range(0, len(self.hfmap_list), 2):
                if i % 2 == 0:
                    _, mask, _ = core.iotools.read_map(self.mask_list[i // 2])
                    mask = utils.set_dim_even(mask)
                uc, arr1, origin = core.iotools.read_map(self.hfmap_list[i])
                uc, arr2, origin = core.iotools.read_map(self.hfmap_list[i + 1])
                arr1 = utils.set_dim_even(arr1) * mask
                arr2 = utils.set_dim_even(arr2) * mask
                if i == 0:
                    nx, ny, nz = arr1.shape
                    map_origin = origin
                    uc_target = uc
                    target_dim = arr1.shape
                    target_pix_size = []
                    for i in range(3):
                        target_pix_size.append(round(uc_target[i] / target_dim[i], 5))
                    #target_pix_size = uc_target[0] / target_dim[0]
                    corner_mask = utils.remove_unwanted_corners(uc, target_dim)
                    if com:
                        com1 = ndimage.measurements.center_of_mass(arr1 * (arr1 >= 0.0))
                        print("COM: ", com1)
                        box_centr = (nx // 2, ny // 2, nz // 2)
                        print(box_centr)
                        self.com1 = com1
                        self.box_centr = box_centr
                        for arr in [arr1, arr2]:
                            arr = shift(arr, np.subtract(box_centr, com1))
                            fhf_lst.append(
                                np.fft.fftshift(
                                    np.fft.fftn(np.fft.fftshift(arr * corner_mask))
                                )
                            )
                    else:
                        for arr in [arr1, arr2]:
                            fhf_lst.append(
                                np.fft.fftshift(
                                    np.fft.fftn(np.fft.fftshift(arr * corner_mask))
                                )
                            )
                else:
                    curnt_pix_size = []
                    for i in range(3):
                        curnt_pix_size.append(round(uc[i] / arr1.shape[i], 5))
                    #curnt_pix_size = uc[0] / arr1.shape[0]
                    for arr in [arr1, arr2]:
                        arr = core.iotools.resample2staticmap(
                            curnt_pix=curnt_pix_size,
                            targt_pix=target_pix_size,
                            targt_dim=target_dim,
                            arr=arr,
                        )
                        if com:
                            com1 = ndimage.measurements.center_of_mass(
                                arr * (arr >= 0.0)
                            )
                            print("COM: ", com1)
                            arr = shift(arr, np.subtract(box_centr, com1))
                        fhf_lst.append(
                            np.fft.fftshift(
                                np.fft.fftn(np.fft.fftshift(arr * corner_mask))
                            )
                        )
            self.map_origin = map_origin
            self.map_unit_cell = uc_target
            self.map_dim = target_dim
            self.fhf_lst = fhf_lst

        if self.mask_list is None:
            for i in range(0, len(self.hfmap_list), 2):
                uc, arr1, origin = core.iotools.read_map(self.hfmap_list[i])
                uc, arr2, origin = core.iotools.read_map(self.hfmap_list[i + 1])
                arr1 = utils.set_dim_even(arr1)
                arr2 = utils.set_dim_even(arr2)
                if i == 0:
                    nx, ny, nz = arr1.shape
                    map_origin = origin
                    uc_target = uc
                    target_dim = arr1.shape
                    target_pix_size = []
                    for i in range(3):
                        target_pix_size.append(round(uc_target[i] / target_dim[i], 5))
                    #target_pix_size = uc_target[0] / target_dim[0]
                    corner_mask = utils.remove_unwanted_corners(uc, target_dim)
                    tmp_lst = []
                    if com:
                        # remove effect from negative values
                        com1 = ndimage.measurements.center_of_mass(arr1 * (arr1 >= 0.0))
                        box_centr = (nx // 2, ny // 2, nz // 2)
                        print("COM, BOX center: ", com1, box_centr)
                        self.com1 = com1
                        self.box_centr = box_centr
                        for arr in [arr1, arr2]:
                            arr = shift(arr, np.subtract(box_centr, com1))
                            tmp_lst.append(arr)
                            fhf_lst.append(
                                np.fft.fftshift(
                                    np.fft.fftn(np.fft.fftshift(arr * corner_mask))
                                )
                            )
                        # reconstruct full map
                        fullmap = ext.half2full.half2full(tmp_lst[0], tmp_lst[1])
                        self.arr_lst.append(fullmap)
                        print(
                            "COM after centering: ",
                            ndimage.measurements.center_of_mass(
                                fullmap * (fullmap >= 0.0)
                            ),
                        )
                        core.iotools.write_mrc(
                            fullmap, "static_centered.mrc", uc_target, map_origin
                        )
                    else:
                        for arr in [arr1, arr2]:
                            tmp_lst.append(arr)
                            fhf_lst.append(
                                np.fft.fftshift(
                                    np.fft.fftn(np.fft.fftshift(arr * corner_mask))
                                )
                            )
                        # reconstruct full map
                        fullmap = ext.half2full.half2full(tmp_lst[0], tmp_lst[1])
                        self.arr_lst.append(fullmap)
                        core.iotools.write_mrc(
                            fullmap, "static_map.mrc", uc_target, map_origin
                        )
                else:
                    tmp_lst = []
                    curnt_pix_size = []
                    for i in range(3):
                        curnt_pix_size.append(round(uc[i] / arr1.shape[i], 5))
                    #curnt_pix_size = uc[0] / arr1.shape[0]
                    for _, arr in enumerate([arr1, arr2]):
                        arr = core.iotools.resample2staticmap(
                            curnt_pix=curnt_pix_size,
                            targt_pix=target_pix_size,
                            targt_dim=target_dim,
                            arr=arr,
                        )
                        if com:
                            com1 = ndimage.measurements.center_of_mass(
                                arr * (arr >= 0.0)
                            )
                            print("COM: ", com1)
                            arr = shift(arr, np.subtract(box_centr, com1))
                            print(
                                "COM after centering: ",
                                ndimage.measurements.center_of_mass(arr * (arr >= 0.0)),
                            )
                        tmp_lst.append(arr)
                        fhf_lst.append(
                            np.fft.fftshift(
                                np.fft.fftn(np.fft.fftshift(arr * corner_mask))
                            )
                        )
                    # reconstruct full map
                    self.arr_lst.append(ext.half2full.half2full(tmp_lst[0], tmp_lst[1]))
            self.map_origin = map_origin
            self.map_unit_cell = uc_target
            self.map_dim = target_dim
            self.fhf_lst = fhf_lst

    def calc_fsc_from_maps(self, fobj):
        # function for only two maps fitting
        nmaps = len(self.fhf_lst)
        fFo_lst = []
        fEo_lst = []
        fBTV_lst = []
        #
        nx, ny, nz = self.fhf_lst[0].shape
        self.nbin, self.res_arr, self.bin_idx = core.restools.get_resolution_array(
            self.map_unit_cell, self.fhf_lst[0]
        )
        for i in range(0, nmaps, 2):
            binfsc, _, _, totalvar, fo, eo = core.fsc.halfmaps_fsc_variance(
                self.fhf_lst[i], self.fhf_lst[i + 1], self.bin_idx, self.nbin
            )
            # weight eo by fsc
            eo = eo * fcodes_fast.read_into_grid(
                self.bin_idx, binfsc, self.nbin, nx, ny, nz
            )
            fFo_lst.append(fo)
            fEo_lst.append(eo)
            fBTV_lst.append(totalvar)
        #
        self.fo_lst = fFo_lst
        self.eo_lst = fEo_lst
        self.totalvar_lst = fBTV_lst


# Likelihood map average
class EmmapAverage:
    def __init__(self, hfmap_list, com=False, phasrand=False, mask_list=None):
        self.hfmap_list = hfmap_list
        self.mask_list = mask_list
        self.map_unit_cell = None
        self.map_origin = None
        self.map_dim = None
        self.com_lst = None
        self.fhf_lst = None
        self.com_lst = None
        self.unmask_fhf_lst = None
        self.phrand_fhf_lst = None
        self.resol_rand = 10.0  # resolution (A) for phase randomisation
        self.com = com
        self.phasrand = phasrand
        self.ceo_lst = None
        self.cbin_idx = None
        self.cdim = None
        self.cbin = None
        self.box_centr = None
        self.s_grid = None
        self.nbin = None
        self.res_arr = None
        self.bin_idx = None
        self.fo_lst = None
        self.eo_lst = None
        self.signalvar_lst = None
        self.totalvar_lst = None
        self.hffsc_lst = None

    def load_maps(self, fobj):
        from scipy import ndimage
        from scipy.ndimage.interpolation import shift
        from numpy.fft import fftn, ifftn, fftshift, ifftshift

        # read masks
        fhf_lst = []
        unmask_fhf_lst = []
        phrand_fhf_lst = []
        com_lst = []
        if self.mask_list is not None:
            if len(self.hfmap_list) // 2 != len(self.mask_list):
                raise SystemExit("mask_list size is not equal to half the size of map_list!")
            for i in range(0, len(self.hfmap_list), 2):
                if i % 2 == 0:
                    _, mask, _ = core.iotools.read_map(self.mask_list[i // 2])
                    mask = utils.set_dim_even(mask)
                uc, arr1, origin = core.iotools.read_map(self.hfmap_list[i])
                uc, arr2, origin = core.iotools.read_map(self.hfmap_list[i + 1])
                arr1 = utils.set_dim_even(arr1)
                arr2 = utils.set_dim_even(arr2)
                if i == 0:
                    fobj.write("Static map:\n")
                    fobj.write(
                        "Mask file: %s\n" % os.path.abspath(self.mask_list[i // 2])
                    )
                    fobj.write("Input files (static map):\n")
                    fobj.write("%s\n" % os.path.abspath(self.hfmap_list[i]))
                    fobj.write("%s\n" % os.path.abspath(self.hfmap_list[i + 1]))
                    nx, ny, nz = arr1.shape
                    assert mask.shape == arr1.shape == arr2.shape
                    if self.com:
                        com1 = ndimage.measurements.center_of_mass(arr1 * mask)
                        com_lst.append(com1)
                        box_centr = (nx // 2, ny // 2, nz // 2)
                        print(box_centr)
                        fobj.write("Center_of_mass coordinates: " + str(com1) + "\n")
                        self.box_centr = box_centr
                        fobj.write(
                            "Center_of_box coordinates: " + str(box_centr) + "\n"
                        )
                    map_origin = origin
                    fobj.write("Map origin coordinates: " + str(map_origin) + "\n")
                    uc_target = uc
                    fobj.write("Unit cell: " + str(uc_target) + "\n")
                    target_dim = arr1.shape
                    fobj.write("Map dimensions: " + str(target_dim) + "\n")
                    target_pix_size = []
                    for i in range(3):
                        target_pix_size.append(round(uc_target[i] / target_dim[i], 5))
                    #target_pix_size = uc_target[0] / target_dim[0]
                    fobj.write("Map pixel size: " + str(target_pix_size) + "\n")
                    corner_mask = utils.remove_unwanted_corners(uc, target_dim)
                    # get resolution grid
                    maxbin = np.amax(np.array([nx // 2, ny // 2, nz // 2]))
                    fobj.write("Creating resolution grid... \n")
                    resol_grid, self.s_grid, _ = fcodes_fast.resolution_grid_full(
                        uc, 0.0, 1, maxbin, nx, ny, nz
                    )
                    if self.phasrand:
                        fobj.write(
                            "Phase randomization using \
                            static half maps. \n"
                        )
                        fobj.write(
                            "Phase randomize resolution:" + str(self.resol_rand) + " \n"
                        )
                    for arr in [arr1, arr2]:
                        arr_mask = arr * mask
                        if self.com:
                            arr_mask = shift(arr_mask, np.subtract(box_centr, com1))
                        fhf_lst.append(fftshift(fftn(fftshift(arr_mask * corner_mask))))
                        print("fhf_lst appended")
                        if self.phasrand:
                            arr_unmask = arr
                            fhf1_randomized = ext.phase_randomize.get_randomized_sf(
                                resol_grid, arr_unmask, self.resol_rand
                            )
                            arr_rand = np.real(ifftn(ifftshift(fhf1_randomized))) * mask
                            if self.com:
                                arr_unmask = shift(
                                    arr_unmask, np.subtract(box_centr, com1)
                                )
                                arr_rand = shift(arr_rand, np.subtract(box_centr, com1))
                            unmask_fhf_lst.append(
                                fftshift(fftn(fftshift(arr_unmask * corner_mask)))
                            )
                            phrand_fhf_lst.append(
                                fftshift(fftn(fftshift(arr_rand * corner_mask)))
                            )

                else:
                    fobj.write("Moving map:\n")
                    fobj.write(
                        "Mask file: %s\n" % os.path.abspath(self.mask_list[i // 2])
                    )
                    fobj.write("Resampling mask...\n")
                    mask = core.iotools.resample2staticmap(
                        target_pix_size, target_dim, uc, mask, fobj=fobj
                    )
                    fobj.write("Input files (moving map):\n")
                    fobj.write("%s\n" % os.path.abspath(self.hfmap_list[i]))
                    fobj.write("%s\n" % os.path.abspath(self.hfmap_list[i + 1]))
                    run_once = True
                    if self.phasrand:
                        fobj.write("Phase randomization using half maps. \n")
                        fobj.write(
                            "Phase randomize resolution:" + str(self.resol_rand) + " \n"
                        )
                    for arr in [arr1, arr2]:
                        fobj.write("Resampling halfmaps...\n")
                        arr = core.iotools.resample2staticmap(
                            target_pix_size, target_dim, uc, arr, fobj=fobj
                        )
                        assert mask.shape == arr.shape
                        arr_mask = arr * mask
                        if self.com:
                            if run_once:
                                com1 = ndimage.measurements.center_of_mass(arr * mask)
                                run_once = False
                            com_lst.append(com1)
                            arr_mask = shift(arr_mask, np.subtract(box_centr, com1))
                        fhf_lst.append(fftshift(fftn(fftshift(arr_mask * corner_mask))))
                        if self.phasrand:
                            arr_unmask = arr
                            fobj.write("Phase randomization...\n")
                            fhf1_randomized = ext.phase_randomize.get_randomized_sf(
                                resol_grid, arr_unmask, self.resol_rand
                            )
                            arr_rand = np.real(ifftn(ifftshift(fhf1_randomized))) * mask
                            if self.com:
                                arr_unmask = shift(
                                    arr_unmask, np.subtract(box_centr, com1)
                                )
                                arr_rand = shift(arr_rand, np.subtract(box_centr, com1))
                            unmask_fhf_lst.append(
                                fftshift(fftn(fftshift(arr_unmask * corner_mask)))
                            )
                            phrand_fhf_lst.append(
                                fftshift(fftn(fftshift(arr_rand * corner_mask)))
                            )
            self.map_origin = map_origin
            self.map_unit_cell = uc_target
            self.map_dim = target_dim
            self.fhf_lst = fhf_lst
            self.com_lst = com_lst
            self.unmask_fhf_lst = unmask_fhf_lst
            self.phrand_fhf_lst = phrand_fhf_lst

        if self.mask_list is None:
            print("Yet to implement")
            exit()
        """if self.mask_list is None: 
            print('Correlation based masks will be generated and used for fitting!')
            from emda import maskmap_class
            obj_maskmap = maskmap_class.MaskedMaps()
            for i in range(0,len(self.hfmap_list),2):
                uc,arr1,origin = read_map(self.hfmap_list[i])
                uc,arr2,origin = read_map(self.hfmap_list[i+1])
                # calculate the mask
                obj_maskmap.generate_mask(arr1, arr2)
                mask = obj_maskmap.mask
                write_mrc(mask,"{0}_{1}.{2}".format('mask',str(i),'mrc'),uc,origin)
                if i == 0:
                    com1 = ndimage.measurements.center_of_mass(arr1 * mask)
                    nx, ny, nz = arr1.shape
                    box_centr = (nx//2, ny//2, nz//2)
                    print(box_centr)
                    self.com1 = com1
                    self.box_centr = box_centr
                    map_origin = origin
                    uc_target = uc
                    target_dim = arr1.shape
                    target_pix_size = uc_target[0]/target_dim[0]
                    # get resolution grid
                    maxbin = np.amax(np.array([nx//2,ny//2,nz//2]))
                    resol_grid, self.s_grid, _ = fcodes_fast.resolution_grid_full(uc,0.0,1,maxbin,nx,ny,nz)
                    #
                    for arr in [arr1, arr2]:
                        arr_unmask = arr
                        fhf1_randomized = get_randomized_sf(resol_grid,arr_unmask,self.resol_rand)
                        arr_rand = np.real(np.fft.ifftn(np.fft.ifftshift(fhf1_randomized))) * mask
                        arr_mask = shift(arr * mask, np.subtract(box_centr,com1))
                        arr_unmask = shift(arr_unmask, np.subtract(box_centr,com1)) 
                        arr_rand = shift(arr_rand, np.subtract(box_centr,com1))
                        fhf_lst.append(np.fft.fftshift(np.fft.fftn(np.fft.fftshift(arr_mask))))
                        unmask_fhf_lst.append(np.fft.fftshift(np.fft.fftn(np.fft.fftshift(arr_unmask))))
                        phrand_fhf_lst.append(np.fft.fftshift(np.fft.fftn(np.fft.fftshift(arr_rand))))
                else:
                    run_once = True
                    for arr in [arr1, arr2]:
                        mask = resample2staticmap(target_pix_size,target_dim,uc,mask)
                        arr = resample2staticmap(target_pix_size,target_dim,uc,arr)
                        arr_unmask = arr
                        fhf1_randomized = get_randomized_sf(resol_grid,arr_unmask,self.resol_rand)
                        arr_rand = np.real(np.fft.ifftn(np.fft.ifftshift(fhf1_randomized))) * mask
                        if run_once: com1 = ndimage.measurements.center_of_mass(arr * mask)
                        run_once = False
                        arr_mask = shift(arr * mask, np.subtract(box_centr,com1))
                        arr_unmask = shift(arr_unmask, np.subtract(box_centr,com1)) 
                        arr_rand = shift(arr_rand, np.subtract(box_centr,com1))
                        fhf_lst.append(np.fft.fftshift(np.fft.fftn(np.fft.fftshift(arr_mask))))
                        unmask_fhf_lst.append(np.fft.fftshift(np.fft.fftn(np.fft.fftshift(arr_unmask))))
                        phrand_fhf_lst.append(np.fft.fftshift(np.fft.fftn(np.fft.fftshift(arr_rand))))
                
            # free memory
            del arr, arr1, arr2
            #            
            self.map_origin     = map_origin
            self.map_unit_cell  = uc_target
            self.map_dim        = target_dim 
            self.fhf_lst        = fhf_lst 
            self.unmask_fhf_lst = unmask_fhf_lst 
            self.phrand_fhf_lst = phrand_fhf_lst"""

    def calc_fsc_variance_from_halfdata(self, fobj):
        nmaps = len(self.fhf_lst)
        fFo_lst = []
        fEo_lst = []
        fBNV_lst = []
        fBSV_lst = []
        fBTV_lst = []
        fBFsc_lst = []
        #
        fobj.write("\n Calculating Fourier Shell Correlation using half data \n")
        fobj.write("   ----------------------------------------------------- \n")
        fobj.write("\n")
        self.nbin, self.res_arr, self.bin_idx = core.restools.get_resolution_array(
            self.map_unit_cell, self.fhf_lst[0]
        )
        idx = np.argmin((self.res_arr - self.resol_rand) ** 2)
        #
        for i in range(0, nmaps, 2):
            # masked data
            fobj.write("\n")
            fobj.write("Masked maps %s\n" % os.path.abspath(self.hfmap_list[i]))
            fobj.write("            %s\n" % os.path.abspath(self.hfmap_list[i + 1]))
            fobj.write("\n")
            (
                bin_fsc,
                noisevar,
                signalvar,
                totalvar,
                fo,
                eo,
            ) = core.fsc.halfmaps_fsc_variance(
                self.fhf_lst[i], self.fhf_lst[i + 1], self.bin_idx, self.nbin
            )
            #
            fobj.write(" bin# \n")
            fobj.write(" bin resolution (A) \n")
            fobj.write(" Signal Variance \n")
            fobj.write(" Noise variance \n")
            fobj.write(" Total variance \n")
            fobj.write(" Halfmap FSC \n")
            fobj.write("\n")
            for j in range(len(self.res_arr)):
                fobj.write(
                    "{:5d} {:6.2f} {:8.4f} {:8.4f} {:8.4f} \
                    {:8.4f}\n".format(
                        j,
                        self.res_arr[j],
                        signalvar[j],
                        noisevar[j],
                        totalvar[j],
                        bin_fsc[j],
                    )
                )
            full_fsc_total = 2.0 * bin_fsc / (1.0 + bin_fsc)
            fFo_lst.append(fo)
            fEo_lst.append(eo)
            fBNV_lst.append(noisevar)
            fBSV_lst.append(signalvar)
            fBTV_lst.append(totalvar)
            fBFsc_lst.append(full_fsc_total)
            if self.phasrand:
                # unmasked data
                fobj.write("Unmasked maps %s\n" % os.path.abspath(self.hfmap_list[i]))
                fobj.write(
                    "              %s\n" % os.path.abspath(self.hfmap_list[i + 1])
                )
                fobj.write("\n")
                (
                    umbin_fsc,
                    umnoisevar,
                    umsignalvar,
                    umtotalvar,
                    _,
                    _,
                ) = core.fsc.halfmaps_fsc_variance(
                    self.unmask_fhf_lst[i],
                    self.unmask_fhf_lst[i + 1],
                    self.bin_idx,
                    self.nbin,
                )
                fobj.write(" bin# \n")
                fobj.write(" bin resolution (A) \n")
                fobj.write(" Signal Variance \n")
                fobj.write(" Noise variance \n")
                fobj.write(" Total variance \n")
                fobj.write(" Halfmap FSC \n")
                fobj.write("\n")
                for j in range(len(self.res_arr)):
                    fobj.write(
                        "{:5d} {:6.2f} {:8.4f} {:8.4f} {:8.4f} \
                        {:8.4f}\n".format(
                            j,
                            self.res_arr[j],
                            umsignalvar[j],
                            umnoisevar[j],
                            umtotalvar[j],
                            umbin_fsc[j],
                        )
                    )
                full_fsc_unmasked = 2.0 * umbin_fsc / (1.0 + umbin_fsc)

                # randomised data
                rbin_fsc, _, _, _, _, _ = core.fsc.halfmaps_fsc_variance(
                    self.phrand_fhf_lst[i],
                    self.phrand_fhf_lst[i + 1],
                    self.bin_idx,
                    self.nbin,
                )
                full_fsc_noise = 2.0 * rbin_fsc / (1.0 + rbin_fsc)
                # fsc_true from Richard's formular
                fsc_true = (full_fsc_total - full_fsc_noise) / (1 - full_fsc_noise)
                # replace fsc_true with fsc_masked_full upto \
                # resol_rand_idx + 2 (RELION uses 2)
                fsc_true[: idx + 2] = full_fsc_total[: idx + 2]

                fobj.write("\n")
                fobj.write(
                    " All halfmap FSCs have been converted \
                    into fullmap FSC \n"
                )
                fobj.write("\n")
                fobj.write(" bin# \n")
                fobj.write(" bin resolution (A) \n")
                fobj.write(" Unmasked FSC \n")
                fobj.write(" Masked FSC \n")
                fobj.write(" Randomized FSC \n")
                fobj.write(" True FSC \n")
                fobj.write("\n")
                for j in range(len(self.res_arr)):
                    fobj.write(
                        "{:5d} {:6.2f} {:8.4f} {:8.4f} {:8.4f} \
                        {:8.4f}\n".format(
                            j,
                            self.res_arr[j],
                            full_fsc_unmasked[j],
                            full_fsc_total[j],
                            full_fsc_noise[j],
                            fsc_true[j],
                        )
                    )

                # plot various FSCs
                fobj.write("\n FSC were plotted into fsc_%d.eps \n" % i)
                core.plotter.plot_nlines(
                    self.res_arr,
                    [full_fsc_unmasked, full_fsc_total, full_fsc_noise, fsc_true],
                    "{0}_{1}.{2}".format("fsc", str(i), "eps"),
                    ["unmasked", "fsc_t", "fsc_n", "fsc_true"],
                )

                fFo_lst.append(fo)
                fEo_lst.append(eo)
                fBNV_lst.append(noisevar)
                fBSV_lst.append(signalvar)
                fBTV_lst.append(totalvar)
                fBFsc_lst.append(fsc_true)
        #
        self.fo_lst = fFo_lst
        self.eo_lst = fEo_lst
        self.signalvar_lst = fBSV_lst
        self.totalvar_lst = fBTV_lst
        self.hffsc_lst = fBFsc_lst

