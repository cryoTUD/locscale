"""
Author: "Rangana Warshamanage, Garib N. Murshudov"
MRC Laboratory of Molecular Biology

This software is released under the
Mozilla Public License, version 2.0; see LICENSE.
"""

# Caller script
from __future__ import absolute_import, division, print_function, unicode_literals
import numpy as np
import fcodes_fast
import os
from emda.config import debug_mode
from emda.core import iotools, maptools, restools, plotter, fsc, quaternions


def read_map(mapname):
    """Reads CCP4 type map (.map) or MRC type map.

    Arguments:
        Inputs:
            mapname: string
                CCP4/MRC map file name

        Outputs:
            uc: float, 1D array
                Unit cell
            arr: float, 3D array
                Map values as Numpy array
            origin: list
                Map origin list
    """
    uc, arr, origin = iotools.read_map(mapname)
    return uc, arr, origin


def read_mtz(mtzfile):
    """Reads mtzfile and returns unit_cell and data in Pandas DataFrame.

    Arguments:
        Inputs:
            mtzfile: string
                MTZ file name

        Outputs:
            uc: float, 1D array
                Unit cell
            df: Pandas data frame
                Map values in Pandas Dataframe
    """
    uc, df = iotools.read_mtz(mtzfile)
    return uc, df


def get_data(struct, resol=5.0, gemmi=False, uc=None, dim=None, maporigin=None):
    """Returns data of a map or a model into an ndarray.

    Reads map data into an ndarray, or if the structure input is an atomic model,
    it calculates the map from the model and returns as an ndarray.

    Arguments:
        Inputs:
            struct: string
                CCP4/MRC map file name or PDB/ENT/CIF file
                resol:  float, optional
                        resolution to calculates map from model. Default is 5.0 A.
                uc: float, 1D array
                    Parameter for modelmap generation. If absent, this will be
                    determined by dim parameter.
                dim: sequence (integers), optional
                    Parameter for modelmap generation. If absent, this will be
                    determined from the size of the molecule.
                maporigin: sequence (integers), optional
                    Parameter for modelmap generation. If present, the calculated map
                    will be shifted according to this information. If absent, this
                    parameter is taken as [0, 0, 0].
                gemmi: bool, optinal
                    Parameter for modelmap calculation. If true the map from coordinates
                    are calculated using Gemmi. Default to False.

        Outputs:
            uc: float, 1D array
                Unit cell
            arr: float, 3D array
                Map values as Numpy array
            origin: list
                Map origin list
    """

    # read map/mrc
    if struct.endswith((".mrc", ".map")):
        uc, arr, orig = read_map(struct)

    # read model pdb/ent/cif
    if struct.endswith((".pdb", ".ent", ".cif")):
        newmodel = struct
        if uc is not None:
            if dim is not None:
                uc = uc
                dim = max(dim)
                orig = [0, 0, 0]
            if dim is None:
                dim = get_dim(model=struct, shiftmodel="new1.cif")
                orig = [0, 0, 0]
                newmodel = "new1.cif"
        if uc is None:
            if dim is not None:
                dim = max(dim)
                uc = np.array([dim, dim, dim, 90.0, 90.0, 90.0], dtype="float")
                orig = [0, 0, 0]
            if dim is None:
                dim = get_dim(model=struct, shiftmodel="new1.cif")
                uc = np.array([dim, dim, dim, 90.0, 90.0, 90.0], dtype="float")
                orig = [-dim // 2, -dim // 2, -dim // 2]  # [0, 0, 0]
                newmodel = "new1.cif"
        if maporigin is None:
            maporigin = orig
        if gemmi:
            modelmap = model2map_gm(
                modelxyz=newmodel,
                dim=[dim, dim, dim],
                resol=resol,
                cell=uc,
                maporigin=maporigin,
            )  
        else:          
            modelmap = model2map(
                modelxyz=newmodel,
                dim=[dim, dim, dim],
                resol=resol,
                cell=uc,
                maporigin=maporigin,
            )
        arr = modelmap
    return uc, arr, orig


def write_mrc(mapdata, filename, unit_cell, map_origin=None, axesorder='ZYX'):
    """Writes 3D Numpy array into MRC file.

    Arguments:
        Inputs:
            mapdata: float, 3D array
                Map values to write
            filename: string
                Output file name
            unit_cell: float, 1D array
                Unit cell params
            map_origin: list, optional
                map origin. Default is [0.0, 0.0, 0.0]
            axesorder: string, optional
                Axes order can be specified for the data to be written.
                By defualt EMDA write data in ZXY convension.
                With this argument, the axes order can be changed.

        Outputs:
            Output MRC file
    """
    iotools.write_mrc(
        mapdata=mapdata, filename=filename, unit_cell=unit_cell, map_origin=map_origin, 
        axesorder=axesorder
    )


def write_mtz(uc, arr, outfile="map2mtz.mtz", resol=None):
    """Writes 3D Numpy array into MTZ file.

    Arguments:
        Inputs:
            uc: float, 1D array
                Unit cell params.
            arr: complex, 3D array
                Map values to write.

        Outputs:
            outfile: string
            Output file name. Default is map2mtz.mtz.
    """
    iotools.write_3d2mtz(unit_cell=uc, mapdata=arr,
                         outfile=outfile, resol=resol)


def resample_data(curnt_pix, targt_pix, arr, targt_dim=None):
    """Resamples a 3D array.

    Arguments:
        Inputs:
            curnt_pix: float list, Current pixel sizes along c, b, a.
            targt_pix: float list, Target pixel sizes along c, b a.
            arr: float, 3D array of map values.
            targt_dim: int list, Target sampling along z, y, x.
            
        Outputs:
            new_arr: float, 3D array
                Resampled 3D array.
    """
    new_arr = iotools.resample2staticmap(
        curnt_pix=curnt_pix, targt_pix=targt_pix, targt_dim=targt_dim, arr=arr
    )
    return new_arr


def estimate_map_resol(hfmap1name, hfmap2name):
    """Estimates map resolution.

    Arguments:
        Inputs:
            hfmap1name: string
                Halfmap 1 name.
            hfmap2name: string
                Halfmap 2 name.

        Outputs:
            map_resol: float
                Map resolution determined by the halfmap FSC.
    """
    map_resol = maptools.estimate_map_resol(
        hfmap1=hfmap1name, hfmap2=hfmap2name)
    return map_resol


def get_map_power(mapname):
    """Calculates the map power spectrum.

    Arguments:
        Inputs:
            mapname: string
                Map name.

        Outputs:
            res_arr: float
                Resolution array.
            power_spectrum: float
                Map power spectrum.
    """
    res_arr, power_spectrum = maptools.get_map_power(mapname)
    return res_arr, power_spectrum


def get_biso_from_model(mmcif_file):
    """Calculates the isotropic B-value of the model.

    Arguments:
        Inputs:
            mmcif_file: string
                mmCIF file name.

        Outputs:
            biso: float
                Model B-iso value.
    """
    biso = maptools.get_biso_from_model(mmcif_file)
    return biso


def get_biso_from_map(halfmap1, halfmap2):
    """Calculates the isotropic B-value of the map.

    Arguments:
        Inputs:
            halfmap1: string
                Halfmap 1 file name.
            halfmap2: string
                Halfmap 2 file name.

        Outputs:
            biso: float
                Map B-iso value.
    """
    biso = maptools.get_biso_from_map(halfmap1=halfmap1, halfmap2=halfmap2)
    return biso


def apply_bfactor_to_map(mapname, bf_arr, mapout):
    """Applies an array of B-factors on the map.

    Arguments:
        Inputs:
            mapname: string
                Map file name.
            bf_arr: float, 1D array
                An array/list of B-factors.
            mapout: bool
                If True, map for each B-factor will be output.

        Outputs:
            all_mapout: complex, ndarray
                4D array containing Fourier coefficients of all maps.
                e.g. all_mapout[:,:,:,i], where i represents map number
                corresponding to the B-factor in bf_arr.
    """
    all_mapout = maptools.apply_bfactor_to_map(
        mapname=mapname, bf_arr=bf_arr, mapout=mapout
    )
    return all_mapout


def map2mtz(mapname, mtzname="map2mtz.mtz", resol=None):
    """Converts a map into MTZ format.

    Arguments:
        Inputs:
            mapname: string
                Map file name.
            mtzname: string
                Output MTZ file name. Default is map2mtz.mtz

        Outputs:
            Outputs MTZ file.
    """
    maptools.map2mtz(mapname=mapname, mtzname=mtzname, resol=resol)


def map2mtzfull(uc, arr1, arr2, mtzname="halfnfull.mtz"):
    """Writes several 3D Numpy arrays into an MTZ file.

    This function accepts densities of two half maps as 3D numpy arrays
    and outputs an MTZ file containing amplitudes of half1, half2 and
    full map. The outfile data labels are H, K, L, Fout1, Fout2, Foutf, Poutf.
    The last four labels correspond to amplitude of halfmap1, amplitudes of
    halfmap2, amplitudes of fullmap and the phase values of fullmap, respectively.

    Arguments:
        Inputs:
            uc: float, 1D array
                Unit cell params.
            arr1: float, 3D array
                Half1 map values.
            arr2: float, 3D array
                Half2 map values.
            mtzname: string, optional
                Output MTZ file name. Default is halfnfull.mtz

        Outputs:
            Outputs an MTZ file containing amplitudes of half maps and
            full map.
    """
    hf1 = np.fft.fftshift(np.fft.fftn(arr1))
    hf2 = np.fft.fftshift(np.fft.fftn(arr2))
    iotools.write_3d2mtz_full(
        unit_cell=uc, hf1data=hf1, hf2data=hf2, outfile=mtzname)


def mtz2map(mtzname, map_size):
    """Converts an MTZ file into MRC format.

    This function converts data in an MTZ file into a 3D Numpy array.
    It combines amplitudes and phases to form Fourier coefficients.

    Arguments:
        Inputs:
            mtzname: string
                MTZ file name.
            map_size: list
                Shape of output 3D Numpy array as a list of three integers.

        Outputs:
            outarr: float
            3D Numpy array of map values.
    """
    arr = maptools.mtz2map(mtzname=mtzname, map_size=map_size)
    return arr


def lowpass_map(uc, arr1, resol, filter="ideal", order=4):
    """Lowpass filters a map to a specified resolution.

    This function applies a lowpass filter on a map to a specified resolution.
    This operations is carried out in the Fourier space. Note that lowpass map
    will have the same shape as input data.

    Arguments:
        Inputs:
            uc: float, 1D array
                Unit cell params.
            arr1: float, 3D array
                3D Numpy array containing map values.
            resol: float
                Resolution cutoff for lowpass filtering in Angstrom units.
            filter: string, optional
                Fiter type to use in truncating Fourier coefficients.
                Currently, only 'ideal' or 'butterworth' filters can be employed.
                Default type is ideal.
            order: integer, optional
                Order of the Butterworth filter. Default is 4.

        Outputs:
            fmap1: complex, 3D array
                Lowpass filtered Fourier coefficeints.
            map1: float, 3D array
                Lowpass filtered map in image/real space
    """
    import emda.ext.lowpass_map as lp

    if filter == "ideal":
        fmap1, map1 = lp.lowpass_map(uc, arr1, resol)
    if filter == "butterworth":
        fmap1, map1 = lp.butterworth(uc, arr1, resol, order)
    return fmap1, map1


def half2full(half1name, half2name, outfile="fullmap.mrc"):
    """Combines half maps to generate full map.

    Arguments:
        Inputs:
            half1name: string
                Name of half map 1.
            half2name: string
                name of half map 2.
            outfile: string, optional
                Name of the output file. Default is fullmap.mrc

        Outputs:
            fullmap: float, 3D array
                3D Numpy array of floats.
    """
    import emda.ext.half2full as h2f

    uc, arr1, origin = read_map(half1name)
    uc, arr2, origin = read_map(half2name)
    fullmap = h2f.half2full(arr1, arr2)
    write_mrc(fullmap, outfile, uc, origin)
    return fullmap


def map_transform(mapname, tra, rot, axr, outname="transformed.mrc"):
    """Imposes a transformation on a map.

    Imposes a transformation (i.e. translation and rotation) on a map
    and returns the transformed map.

    Arguments:
        Inputs:
            mapname: string
                Name of the input map.
            tra: list of three floats values
                Translation vector as a list in Angstrom units.
            rot: float
                Rotation to apply in degrees.
            axr: list of three integers
                Axis to rotation. e.g [1, 0, 0]
            outname: string, optional
                Name of the transformed map. Default is transformed.mrc.

        Outputs:
            transformed_map: float, 3D array
                3D Numpy array of floats.
    """
    from emda.ext import transform_map

    transformed_map = transform_map.map_transform(
        mapname, tra, rot, tuple(axr), outname
    )
    return transformed_map


def halfmap_fsc(half1name, half2name, filename=None, maskname=None):
    """Computes Fourier Shell Correlation (FSC) using half maps.

    Computes Fourier Shell Correlation (FSC) using half maps.
    FSC is not corrected for mask effect in this implementation.

    Arguments:
        Inputs:
            half1name: string
                Name of the half map 1.
            half2name: string
                Name of the half map 2.
            filename: string
                If present, statistics will be printed into this file.
            maskname: String
                If present, input maps will be masked before computing FSC.

        Outputs:
            res_arr: float, 1D array
                Linear array of resolution in Angstrom units.
            bin_fsc: float, 1D array
                Linear array of FSC in each resolution bin.
    """
    import os

    uc, arr1, _ = iotools.read_map(half1name)
    uc, arr2, _ = iotools.read_map(half2name)
    hf1 = np.fft.fftshift(np.fft.fftn(arr1))
    hf2 = np.fft.fftshift(np.fft.fftn(arr2))
    fsc_list = []
    nbin, res_arr, bin_idx = restools.get_resolution_array(uc, hf1)
    (
        _,
        _,
        noisevar,
        signalvar,
        totalvar,
        bin_fsc,
        bincount,
    ) = fcodes_fast.calc_fsc_using_halfmaps(
        hf1, hf2, bin_idx, nbin, debug_mode, hf1.shape[0], hf1.shape[1], hf1.shape[2]
    )
    fsc_list.append(bin_fsc)
    if maskname is not None:
        _, mask, _ = read_map(maskname)
        arr1_msk = arr1 * mask
        arr2_msk = arr2 * mask
        hf1msk = np.fft.fftshift(np.fft.fftn(arr1_msk))
        hf2msk = np.fft.fftshift(np.fft.fftn(arr2_msk))
        (
            _,
            _,
            msk_noisevar,
            msk_signalvar,
            msk_totalvar,
            msk_bin_fsc,
            msk_bincount,
        ) = fcodes_fast.calc_fsc_using_halfmaps(
            hf1msk,
            hf2msk,
            bin_idx,
            nbin,
            debug_mode,
            hf1msk.shape[0],
            hf1msk.shape[1],
            hf1msk.shape[2],
        )
        fsc_list.append(msk_bin_fsc)
    if filename is not None:
        tdata = open(filename, "w")
        tdata.write("halfmap1 file: %s\n" % os.path.abspath(half1name))
        tdata.write("halfmap2 file: %s\n" % os.path.abspath(half2name))
        tdata.write("\n")
        tdata.write("***** Unmasked statistics *****\n")
        tdata.write("\n")
        tdata.write("bin # \n")
        tdata.write("resolution (Ang.) \n")
        tdata.write("signal variance \n")
        tdata.write("noise variance \n")
        tdata.write("total variance \n")
        tdata.write("halfmap fsc \n")
        tdata.write("# reflx \n")
        i = -1
        for sv, nv, tv, fsci, nfc in zip(
            signalvar, noisevar, totalvar, bin_fsc, bincount
        ):
            i += 1
            tdata.write(
                "{:-3d} {:-6.2f} {:-14.4f} {:-14.4f} {:-14.4f} {:-14.4f} {:-10d}\n".format(
                    i, res_arr[i], sv, nv, tv, fsci, nfc
                )
            )
        if maskname is not None:
            tdata.write("\n")
            tdata.write("***** Masked statistics *****\n")
            tdata.write("\n")
            tdata.write("bin # \n")
            tdata.write("resolution (Ang.) \n")
            tdata.write("signal variance \n")
            tdata.write("noise variance \n")
            tdata.write("total variance \n")
            tdata.write("halfmap fsc \n")
            tdata.write("# reflx \n")
            i = -1
            for sv, nv, tv, fsci, nfc in zip(
                msk_signalvar, msk_noisevar, msk_totalvar, msk_bin_fsc, msk_bincount
            ):
                i += 1
                tdata.write(
                    "{:-3d} {:-6.2f} {:-14.4f} {:-14.4f} {:-14.4f} {:-14.4f} {:-10d}\n".format(
                        i, res_arr[i], sv, nv, tv, fsci, nfc
                    )
                )
    if maskname is not None:
        print("Bin    Resolution     unmasked-FSC   masked-FSC")
        for i in range(len(res_arr)):
            print(
                "{:5d} {:6.2f} {:14.4f} {:14.4f}".format(
                    i, res_arr[i], bin_fsc[i], msk_bin_fsc[i]
                )
            )
    else:
        print("Bin    Resolution     FSC")
        for i in range(len(res_arr)):
            print("{:5d} {:6.2f} {:14.4f}".format(i, res_arr[i], bin_fsc[i]))
    return res_arr, fsc_list


def halfmap_fsc_ph(half1name, half2name, filename="halffsc.txt", maskname=None):
    """Computes Fourier Shell Correlation (FSC) using half maps.

    Computes Fourier Shell Correlation (FSC) using half maps.
    FSC is not corrected for mask effect in this implementation.

    Arguments:
        Inputs:
            half1name: string
                Name of the half map 1.
            half2name: string
                Name of the half map 2.
            filename: string
                If present, statistics will be printed into this file.
            maskname: String
                If present, input maps will be masked before computing FSC.

        Outputs:
            res_arr: float, 1D array
                Linear array of resolution in Angstrom units.
            bin_fsc: float, 1D array
                Linear array of FSC in each resolution bin.
    """
    import os

    tdata = open(filename, "w")

    uc, arr1, _ = iotools.read_map(half1name)
    uc, arr2, _ = iotools.read_map(half2name)
    hf1 = np.fft.fftshift(np.fft.fftn(arr1))
    hf2 = np.fft.fftshift(np.fft.fftn(arr2))
    fsc_list = []
    nbin, res_arr, bin_idx = restools.get_resolution_array(uc, hf1)
    (
        _,
        _,
        noisevar,
        signalvar,
        totalvar,
        bin_fsc,
        bincount,
    ) = fcodes_fast.calc_fsc_using_halfmaps(
        hf1, hf2, bin_idx, nbin, debug_mode, hf1.shape[0], hf1.shape[1], hf1.shape[2]
    )
    fsc_list.append(bin_fsc)
    tdata.write("halfmap1 file: %s\n" % os.path.abspath(half1name))
    tdata.write("halfmap2 file: %s\n" % os.path.abspath(half2name))
    tdata.write("\n")
    tdata.write("***** Unmasked statistics *****\n")
    tdata.write("\n")
    tdata.write("bin # \n")
    tdata.write("resolution (Ang.) \n")
    tdata.write("signal variance \n")
    tdata.write("noise variance \n")
    tdata.write("total variance \n")
    tdata.write("halfmap fsc \n")
    tdata.write("# reflx \n")
    i = -1
    for sv, nv, tv, fsci, nfc in zip(signalvar, noisevar, totalvar, bin_fsc, bincount):
        i += 1
        tdata.write(
            "{:-3d} {:-6.2f} {:-14.4f} {:-14.4f} {:-14.4f} {:-14.4f} {:-10d}\n".format(
                i, res_arr[i], sv, nv, tv, fsci, nfc
            )
        )
    if maskname is not None:
        _, mask, _ = read_map(maskname)
        from emda.ext.phase_randomize import phase_randomized_fsc

        idx = np.argmin((bin_fsc - 0.8) ** 2)
        resol_rand = res_arr[idx]
        fsc_list_ph, msk_bincount = phase_randomized_fsc(
            arr1=arr1,
            arr2=arr2,
            mask=mask,
            bin_idx=bin_idx,
            res_arr=res_arr,
            fobj=tdata,
            # resol_rand=resol_rand,
        )
        fsc_list.append(fsc_list_ph[0])
        fsc_list.append(fsc_list_ph[1])
        fsc_list.append(fsc_list_ph[2])

        tdata.write("\n")
        tdata.write("***** Fourier Shell Correlation *****\n")
        tdata.write("\n")
        tdata.write("bin # \n")
        tdata.write("resolution (Ang.) \n")
        tdata.write("Unmask FSC \n")
        tdata.write("Masked FSC \n")
        tdata.write("Noise FSC \n")
        tdata.write("True FSC \n")
        tdata.write("# reflx \n")
        i = -1
        for umf, mf, nf, tf, nfc in zip(
            fsc_list[0], fsc_list[1], fsc_list[2], fsc_list[3], msk_bincount
        ):
            i += 1
            tdata.write(
                "{:-3d} {:-6.2f} {:-14.4f} {:-14.4f} {:-14.4f} {:-14.4f} {:-10d}\n".format(
                    i, res_arr[i], umf, mf, nf, tf, nfc
                )
            )
    if len(fsc_list) == 4:
        plotter.plot_nlines(
            res_arr,
            fsc_list,
            "halfmap_fsc_ph.eps",
            curve_label=["Unmask", "Masked", "Noise", "Corrected"],
            plot_title="Halfmap FSC",
        )
    return res_arr, fsc_list


def get_variance(half1name, half2name, filename=None, maskname=None):
    """Returns noise and signal variances of half maps.

    Returns noise and signal variances of half maps. Return values are not
    corrected for full map.

    Arguments:
        Inputs:
            half1name: string
                Name of the half map 1.
            half2name: string
                Name of the half map 2.
            filename: string
                If present, statistics will be printed into this file.
            maskname: String
                If present, input maps will be masked before computing variances.

        Outputs:
            res_arr: float, 1D array
                Linear array of resolution in Angstrom units.
            noisevar: float, 1D array
                Linear array of noise variance in each resolution bin.
            signalvar: float, 1D array
                Linear array of signal variance in each resolution bin.
    """
    import os

    uc, arr1, _ = iotools.read_map(half1name)
    uc, arr2, _ = iotools.read_map(half2name)
    if maskname is not None:
        _, mask, _ = read_map(maskname)
        arr1 = arr1 * mask
        arr2 = arr2 * mask
    hf1 = np.fft.fftshift(np.fft.fftn(arr1))
    hf2 = np.fft.fftshift(np.fft.fftn(arr2))
    nbin, res_arr, bin_idx = restools.get_resolution_array(uc, hf1)
    (
        _,
        _,
        noisevar,
        signalvar,
        totalvar,
        bin_fsc,
        bincount,
    ) = fcodes_fast.calc_fsc_using_halfmaps(
        hf1, hf2, bin_idx, nbin, debug_mode, hf1.shape[0], hf1.shape[1], hf1.shape[2]
    )
    if filename is not None:
        tdata = open(filename, "w")
        tdata.write("halfmap1 file: %s\n" % os.path.abspath(half1name))
        tdata.write("halfmap2 file: %s\n" % os.path.abspath(half2name))
        tdata.write("\n")
        tdata.write("bin # \n")
        tdata.write("resolution (Ang.) \n")
        tdata.write("signal variance \n")
        tdata.write("noise variance \n")
        tdata.write("total variance \n")
        tdata.write("halfmap fsc \n")
        tdata.write("# reflx \n")
        i = -1
        for sv, nv, tv, fsci, nfc in zip(
            signalvar, noisevar, totalvar, bin_fsc, bincount
        ):
            i += 1
            tdata.write(
                "{:-3d} {:-6.2f} {:-14.4f} {:-14.4f} {:-14.4f} {:-14.4f} {:-10d}\n".format(
                    i, res_arr[i], sv, nv, tv, fsci, nfc
                )
            )
    return res_arr, noisevar, signalvar


def twomap_fsc(map1name, map2name, fobj=None, xmlobj=None):
    """Returns Fourier Shell Correlation (FSC) between any two maps.

    Computes Fourier Shell Correlation (FSC) using any two maps.

    Arguments:
        Inputs:
            map1name: string
                Name of the map 1.
            map2name: string
                Name of the map 2.
            fobj: file object for logging
                If present, statistics will be printed into this file.
            xmlobj: xml object
                If present, statistics will be printed into an XML file.

        Outputs:
            res_arr: float, 1D array
                Linear array of resolution in Angstrom units.
            bin_fsc: float, 1D array
                Linear array of FSC in each resolution bin.
    """
    import os

    uc, arr1, _ = iotools.read_map(map1name)
    uc, arr2, _ = iotools.read_map(map2name)
    f1 = np.fft.fftshift(np.fft.fftn(arr1))
    f2 = np.fft.fftshift(np.fft.fftn(arr2))
    nbin, res_arr, bin_idx = restools.get_resolution_array(uc, f1)
    bin_fsc = fsc.anytwomaps_fsc_covariance(f1, f2, bin_idx, nbin)[0]
    if xmlobj is not None:
        xmlobj.map1path = os.path.abspath(map1name)
        xmlobj.map2path = os.path.abspath(map2name)
        xmlobj.res_arr = res_arr
        xmlobj.fsc = bin_fsc
        # xmlobj.outmap = "fullmap.mtz" # this file may not present
        xmlobj.write_xml()
    if fobj is not None:
        fobj.write("map1 file: %s\n" % os.path.abspath(map1name))
        fobj.write("map2 file: %s\n" % os.path.abspath(map2name))
        fobj.write("\n")
        fobj.write("bin # \n")
        fobj.write("resolution (Ang.) \n")
        fobj.write("fsc \n")
        for ibin, fsc1 in enumerate(bin_fsc):
            fobj.write("{:5d} {:6.2f} {:6.3f}\n".format(
                ibin, res_arr[ibin], fsc1))
    print("Bin      Resolution     FSC")
    for ibin, fsc2 in enumerate(bin_fsc):
        print("{:5d} {:6.2f} {:6.3f}".format(ibin, res_arr[ibin], fsc2))
    return res_arr, bin_fsc


def balbes_data(map1name, map2name, fsccutoff=0.5, mode="half"):
    """Returns data required for Balbes pipeline.

    Required data is output with their references in EMDA.xml.

    Arguments:
        Inputs:
            map1name: string
                Name of the map 1.
            map2name: string
                Name of the map 2.
            fsccutoff: float, optional
                FSC of desired resolution. Defualt is 0.5
            mode: string
                Mode can be either 'half' or 'any'. If the input maps are
                the half maps, mode should be 'half'. Otherwise, mode should be 'any'.
                Default mode is half.

        Outputs:
            Outputs EMDA.xml containing data and references to other data.
            No return variables.
    """
    from emda.ext import xmlclass

    if mode == "half":
        prepare_refmac_data(
            hf1name=map1name, hf2name=map2name, fsccutoff=fsccutoff)
    else:
        xml = xmlclass.Xml()
        res_arr, bin_fsc = twomap_fsc(
            map1name=map1name, map2name=map2name, xmlobj=xml)


def singlemap_fsc(map1name, knl=3):
    """Returns Fourier Shell Correlation (FSC) of a map.

    Computes Fourier Shell Correlation (FSC) between a map and its
    reconstituted other half from neighbough Fourier coefficients.
    This method can be used to estimate FSC based resolution. However,
    results seem to be reliable when an unfiltered map is used.

    Arguments:
        Inputs:
            map1name: string
                Name of the map.
            knl: integer, optional
                Radius of the integrating kernel.

        Outputs:
            res_arr: float, 1D array
                Linear array of resolution in Angstrom units.
            bin_fsc: float, 1D array
                Linear array of FSC in each resolution bin.
            Outputs reconstituted map as 'fakehalf.mrc'
    """
    from scipy import interpolate
    from emda.ext import fakehalf
    from emda.ext.mapfit.mapaverage import set_array

    uc, arr1, origin = iotools.read_map(map1name)
    f1 = np.fft.fftshift(np.fft.fftn(np.fft.fftshift(arr1)))
    f2 = fakehalf.fakehalf(f_map=f1, knl=knl)
    data2write = np.real(np.fft.ifftshift(np.fft.ifftn(np.fft.ifftshift(f2))))
    write_mrc(data2write, "fakehalf.mrc", uc, origin)
    nbin, res_arr, bin_idx = restools.get_resolution_array(uc, f1)
    bin_fsc, _, _, _, _, _ = fsc.halfmaps_fsc_variance(f1, f2, bin_idx, nbin)
    print("Resolution bin     FSC")
    for i, _ in enumerate(res_arr):
        print("{:.2f} {:.4f}".format(res_arr[i], bin_fsc[i]))
    # deciding resolution
    bin_fsc_trunc = set_array(bin_fsc, thresh=0.15)
    dist05 = np.sqrt((bin_fsc_trunc - 0.5) ** 2)
    indx = fakehalf.get_index(dist05)
    lim1 = bin_fsc[indx]
    res1 = res_arr[indx]
    if lim1 > 0.5:
        lim2 = bin_fsc[indx + 1]
        res2 = res_arr[indx + 1]
        fsc_seq = [lim1, lim2]
        res_seq = [res1, res2]
    elif lim1 < 0.5:
        lim2 = bin_fsc[indx - 1]
        res2 = res_arr[indx - 1]
        fsc_seq = [lim2, lim1]
        res_seq = [res2, res1]
    f = interpolate.interp1d(fsc_seq, res_seq)
    map_resol = f(0.5)
    print("Map resolution (A): ", map_resol)
    return res_arr, bin_fsc, map_resol


def get_fsc(arr1, arr2, uc):
    """Returns FSC as a function of resolution

    Arguments:
        Inputs:
            arr1: float, ndarray
                Density array 1.
            arr2: float, ndarray
                Density array 2.
            uc: float, 1D array
                Unit cell

        Outputs:
            res_arr: float, 1D array
                Linear array of resolution in Angstrom units.
            bin_fsc: float, 1D array
                Linear array of FSC in each resolution bin.
    """
    fmap1 = np.fft.fftshift(np.fft.fftn(arr1))
    fmap2 = np.fft.fftshift(np.fft.fftn(arr2))
    nbin, res_arr, bin_idx = restools.get_resolution_array(uc, fmap1)
    bin_fsc = fsc.anytwomaps_fsc_covariance(fmap1, fmap2, bin_idx, nbin)[0]
    return res_arr, bin_fsc


def mask_from_halfmaps(uc, half1, half2, radius=4, iter=1, dthresh=None):
    """Generates a mask from half maps.

    Generates a mask from half maps based on real space local correlation.

    Arguments:
        Inputs:
            uc: float, 1D array
                Unit cell parameters.
            half1: float, 3D array
                Half map 1 data.
            half2: float, 3D array
                Half map 2 data.
            radius: integer, optional
                Radius of integrating kernel in voxels. Default is 4.
            iter: integer,optional
                Number of dilation cycles. Default is 1 cycle.
            dthresh: float, optional
                The halfmap densities will be thresholded at this value prior
                calculating local correlation. If the value is not given,
                EMDA takes this values as the value at which the cumulative
                density probability is 0.99. This is the default. However,
                it is recommomded that this value should be supplied
                by the user and proven to be useful for cases those have
                micellular densities around protein.
                

        Outputs:
            mask: float, 3D array
                3D Numpy array of correlation mask.
    """
    from emda.ext import realsp_local, maskmap_class

    arr1, arr2 = half1, half2
    # normalized maps
    #if norm:
    #    hf1 = np.fft.fftshift(np.fft.fftn(arr1))
    #    hf2 = np.fft.fftshift(np.fft.fftn(arr2))
    #    arr1, arr2 = realsp_local.hfdata_normalized(hf1=hf1, hf2=hf2, uc=uc)
    obj_maskmap = maskmap_class.MaskedMaps()
    obj_maskmap.smax = radius
    obj_maskmap.arr1 = arr1
    obj_maskmap.arr2 = arr2
    obj_maskmap.uc = uc
    obj_maskmap.iter = iter
    obj_maskmap.dthresh = dthresh
    obj_maskmap.generate_mask()
    return obj_maskmap.mask


def mask_from_map(
    uc,
    arr,
    kern=5,
    resol=15,
    filter="butterworth",
    order=1,
    prob=0.99,
    itr=3,
    orig=None,
):
    """Generates a mask from a map.

    Generates a mask from a map.

    Arguments:
        Inputs:
            uc: float, 1D array
                Unit cell parameters.
            arr: float, 3D array
                Map data.
            half2: float, 3D array
                Half map 2 data.
            kern: integer, optional
                Radius of integrating kernel in voxels. Default is 5.
            resol: float, optional
                Resolution cutoff for lowpass filtering in Angstrom units.
                Default is 15 Angstrom.
            filter: string,optional
                Filter type to use with lowpass filtering. Default is butterworth.
            order: integer, optional
                Butterworth filter order. Default is 1.
            prob: float, optional
                Cumulative probability cutoff to decide the density threshold.
                Default value is 0.99.
            itr: integer, optional
                Number of dilation cycles. Default is 3 cycles.
            orig: list of three integer values.
                Map origin. e.g. [0, 0, 0]

        Outputs:
            mask: float, 3D array
                3D Numpy array of the mask.
            Outputs lowpass.mrc and mapmask.mrc files.
    """
    from emda.ext import maskmap_class

    _, arrlp = lowpass_map(uc, arr, resol, filter, order=order)
    # write_mrc(arrlp, "lowpass.mrc", uc, orig)
    mask = maskmap_class.mapmask(
        arr=arrlp, uc=uc, kern_rad=kern, prob=prob, itr=itr)
    write_mrc(mask, "mapmask.mrc", uc, orig)
    return mask


def mask_from_atomic_model(mapname, modelname, atmrad=5, binary_mask=False):
    """Generates a mask from atomic coordinates.

    Generates a mask from coordinates. First, atomic positions are
    mapped onto a 3D grid. Second, each atomic position is convluted
    with a sphere whose radius is defined by the atmrad paramter.
    Next, one pixel layer dialtion followed by the smoothening of
    edges.

    Arguments:
        Inputs:
            mapname: string
                Name of the map file. This is needed to get the
                sampling, unit cell and origin for the new mask.
                Allowed formats are - MRC/MAP
            modelname: string
                Atomic model name. Allowed formats are - PDB/CIF
            atmrad: float
                Radius of the sphere to be placed on atomic positions in Angstroms.
                Default is 5 A.
            binary_mask: bool
                If used, binary mask will also be output. Default to False

        Outputs:
            mask: float, 3D array
                3D Numpy array of the mask.
            Outputs emda_model_mask.mrc.
    """
    from emda.ext.maskmap_class import mask_from_coordinates

    return mask_from_coordinates(mapname=mapname, 
                                 modelname=modelname, 
                                 atmrad=atmrad, 
                                 binary_mask=binary_mask)


def sphere_kernel_softedge(radius=5):
    """Generates a soft-edged spherical kernel.

    Arguments:
        Inputs:
            radius: integer, optional
                Radius of integrating kernel in voxels. Default is 5.

        Outputs:
            kernel: float, 3D array
                3D Numpy array of spherical kernel.
    """
    kernel = restools.create_soft_edged_kernel_pxl(radius)
    return kernel


def overlay_maps(
    maplist,
    rot=[0.0],
    ncy=100,
    modelres=None,
    masklist=None,
    modellist=None,
    tra=None,
    axr=None,
    fitres=None,
    nocom=False,
):
    """Superimposes several maps.

    Superimposes several maps using a likelihood function. All maps are
    overlaid on the first map.

    Arguments:
        Inputs:
            maplist: list of strings
                List of maps to overlay. The first of the list
                is take as static and it has to be an EM map. The rest can either
                be maps or coordinate files [to combine with modelres argument].
                Coordinate based maps are calculated using REFMAC for fitting. 
            masklist: list of strings
                List of masks to apply on maps.
            rot: float list, optional
                List of initial rotations in degrees. Default to 0.0.
                This list must not include the rotation for the first map to itself
                (because it is zero). Rotations of each map relative to the first
                map should be give, eg. [1., 2., 5.]. If only one value is give [1.],
                that value will be copied to form the full set eg. [1., 1., 1.]
            axr: list of integers, optional
                list of rotation axes. Default to [1, 0, 0].
                Rotation axis for each map except the static map (first map) should eb give.
                If only one axis is specified, that will be copied to form the full list.
            tra: list, optional
                List of translation vector(z, y, x) in Angstroms. Default to [0.0, 0.0, 0.0]
                If only one vector is give, that will be copied to form the full list.
            ncy: integer, optional
                Number of fitting cycles. Default to 100.
            fitres: float, optional
                Maximum resolution for map overlay. Default to 0.0 A (i.e. all data)
            modelres: float, optional
                Resolution for model-based map calculation if there are 
                coordinate files included in the maplist.
                In such case, this resolution is mandatory. All maps will be calculated to
                this resolution.
            modellist: list of strings, optional
                List of coornate files to apply the transformation estimated
                by map overlay. 

        Outputs:
            Outputs a series of overlaid maps (fitted_map_x.mrc).
    """
    #from emda.ext.mapfit import mapoverlay
    from emda.ext.overlay import overlay

    if axr is None:
        axr = [1, 0, 0]
    if tra is None:
        tra = [0., 0., 0.]
    axes = [axr[x:x+3] for x in range(0, len(axr), 3)]
    translations = [tra[x:x+3] for x in range(0, len(tra), 3)]
    rotations = list(rot)
    if len(axes) == 1:
        # duplicate axes
        l1 = [x for _ in range(len(maplist)-1) for x in axes[0]]
        axes = [l1[x:x+3] for x in range(0, len(l1), 3)]
    if len(translations) == 1:
        # duplicate translations
        t1 = [x for _ in range(len(maplist)-1) for x in translations[0]]
        translations = [t1[x:x+3] for x in range(0, len(t1), 3)]
    if len(rotations) == 1:
        rotations = np.zeros(len(maplist)-1, 'float') + rotations[0]
    assert len(maplist) - 1 == len(axes) == len(rotations) == len(translations)
    qlist = []
    for axis, ang in zip(axes, rotations):
        qlist.append(quaternions.get_quaternion([tuple(axis), ang]))
    fobj = open("EMDA_overlay.txt", "w")
    emmap1, rotmat_list, trans_list = overlay(
        maplist=maplist,
        ncycles=ncy,
        tlist=translations,
        qlist=qlist,
        masklist=masklist,
        modellist=modellist,
        fobj=fobj,
        fitres=fitres,
        modelres=modelres,
        nocom=nocom,
    )
    return emmap1, rotmat_list, trans_list

def average_maps(
    maplist,
    rot=0.0,
    ncy=5,
    res=6,
    interp="linear",
    fit=True,
    tra=None,
    axr=None,
    fobj=None,
    masklist=None,
):
    """Calculates the best average maps using Bayesian principles.

    Calculates the best average map using Bayesian principles. This is done in two steps;
    1. Parameter estimation using a likelihood function, 2. Best map calculation.
    Parameter estimation is similar to map overlay where each map is brought onto
    static map by maximizing the overlap. The best maps are calculated using
    superimposed maps.

    Arguments:
        Inputs:
            maplist: list
                List of half maps to average.
            masklist: list, optional
                List of masks to apply on maps. len(masklist) == len(maplist) // 2
            rot: float, optional
                Initial rotation in degrees. Default is 0.0.
            axr: list, optional
                Rotation axis. Default is [1, 0, 0].
            tra: list, optional
                Translation vector in fractional units. Default is [0.0, 0.0, 0.0]
            res: float, optional
                Fit start resolution in Angstrom units. Default is 6.0 Angstrom.
            ncy: integer, optional
                Number of fitting cycles. Default is 5.
            interp: string, optional
                Interpolation type either "linear" or "cubic".
                Default is linear.
            fobj: string
                File object for logging. If None given, EMDA_average.txt will be output.
            fit: bool, optional
                If True, map fitting will be carried out before average map calculation.
                Default is True.

        Outputs:
            Outputs a series of overlaid maps (fitted_map_?.mrc).
            Also, outputs a series of average maps (avgmap_?.mrc)
    """
    from emda.ext.mapfit import mapaverage

    if axr is None:
        axr = [1, 0, 0]
    if tra is None:
        tra = [0.0, 0.0, 0.0]
    if fobj is None:
        fobj = open("EMDA_average.txt", "w")
    theta_init = [tuple(axr), rot]
    mapaverage.main(
        maplist=maplist,
        masklist=masklist,
        ncycles=ncy,
        t_init=tra,
        theta_init=theta_init,
        smax=res,
        fobj=fobj,
        interp=interp,
        fit=True,
    )


def realsp_correlation(
    half1map,
    half2map,
    kernel_size=5,
    norm=False,
    model=None,
    model_resol=None,
    mask_map=None,
    lgf=None,
):
    """Calculates local correlation in real/image space.

    Arguments:
        Inputs:
            half1map: string
                Name of half map 1.
            half1map: string
                Name of half map 2.
            kernel_size: integer, optional
                Radius of integration kernal in pixels. Default is 5.
            norm: bool, optional
                If True, correlation will be carried out on normalized maps.
                Default is False.
            model: string, optional
                An argument for model based map calculation using REFMAC.
                Name of model file (cif/pdb). If present, map-model local
                correlation will be calculated.
            model_resol: float, optional
                An argument for model based map calculation using REFMAC.
                Resolution to calculate model based map. If absent, FSC based
                resolution cutoff will be employed.
            mask_map: string, optional
                Mask file to apply on correlation maps. If not given, correlation based
                mask will be employed.
            lgf: string, optional
                An argument for model based map calculation using REFMAC.
                Ligand description file (cif).

        Outputs:
            Following maps are written out:
            rcc_halfmap_smax?.mrc - reals space half map local correlation.
            rcc_fullmap_smax?.mrc - correlation map corrected to full map
                using the formula 2 x FSC(half) / (1 + FSC(half)).
            If a model included, then
            rcc_mapmodel_smax?.mrc - local correlation map between model and
                full map.
            rcc_truemapmodel_smax?.mrc - truemap-model correaltion map for
                validation purpose.
    """
    from emda.ext import realsp_local

    rcc = realsp_local.RealspaceLocalCC()
    rcc.hfmap1name = half1map
    rcc.hfmap2name = half2map
    rcc.kern_rad = kernel_size
    rcc.model = model
    rcc.model_resol = model_resol
    rcc.maskname = mask_map
    rcc.lgf = lgf
    rcc.norm = norm
    rcc.rcc()


def b_from_correlation(
    half1map,
    half2map,
    resol,
    kernel_size=5,
    mask_map=None,
):
    from emda.ext import realsp_local

    realsp_local.bfromcc(
        half1_map=half1map,
        half2_map=half2map,
        kernel_size=kernel_size,
        resol=resol,
        mask_map=mask_map,
    )


def realsp_correlation_mapmodel(
    fullmap,
    model,
    resol,
    kernel_size=5,
    norm=False,
    mask_map=None,
    lgf=None,
):
    """Calculates real space local correlation between map and model.

    Arguments:
        Inputs:
            fullmap: string
                Name of the map.
            model: string
                An argument for model based map calculation using REFMAC.
                Name of model file (cif/pdb/ent/mtz/mrc).
            resol: float
                An argument for model based map calculation using REFMAC.
                Resolution to calculate model based map.
                All maps are lowpass filtered to this resolution
                before calculating local correlation.
            kernel_size: integer, optional
                Radius of integration kernal in pixels. Default is 5.
            mask_map: string, optional
                Mask file to apply on correlation maps.
            norm: bool, optional
                Defalt to False.
                If True, correlation will be carried out on normalized maps.
                Default is False.
            lgf: string, optional
                An argument for model based map calculation using REFMAC.
                Ligand description file (cif).

        Outputs:
            Following maps are written out:
            modelmap.mrc - model based map (only if atomic model is given).
            mmcc_mapmodel.mrc - real space local correlation map.
    """
    from emda.ext import realsp_local

    realsp_local.mapmodel_rcc(
        fullmap=fullmap,
        model=model,
        kernel_size=kernel_size,
        resol=resol,
        mask_map=mask_map,
        norm=norm,
        lgf=lgf,
    )


def fouriersp_correlation(half1_map, half2_map, kernel_size=5, mask=None):
    """Calculates Fourier space local correlation using half maps.

    Arguments:
        Inputs:
            half1_map: string
                Name of half map 1.
            half2_map: string
                Name of half map 2.
            kernel_size: integer, optional
                Radius of integration kernal. Default is 5.

        Outputs:
            Following maps are written out:
            fouriercorr3d_halfmaps.mrc - local correlation in half maps.
            fouriercorr3d_fullmap.mrc - local correlation in full map
                using the formula 2 x FSC(half) / (1 + FSC(half)).
            fouriercorr3d_truemap.mrc - local correlation in true map.
                Useful for validation purpose.
    """
    from emda.ext import fouriersp_local

    fouriersp_local.fcc(
        half1_map=half1_map, half2_map=half2_map, kernel_size=kernel_size, maskmap=mask
    )


def map_model_validate(
    half1map,
    half2map,
    modelfpdb,
    bfac=0.0,
    lig=True,
    model1pdb=None,
    mask=None,
    modelresol=None,
    lgf=None,
):
    """Calculates various FSCs for maps and model validation.

    Arguments:
        Inputs:
            half1map: string
                Name of half map 1.
            half2map: string
                Name of half map 2.
            modelfpdb: string
                Name of the model refined against full map in cif/pdb/ent
                formats.
            model1pdb: string, optional
                Name of the model refined against one half map in cif/pdb/ent
                formats. If included, FSC between that and half maps will be
                calculated.
            mask: string, optional
                Name of the mask file. It will apply on half maps before
                computing FSC. If not included, a correlation based masked will
                employed.
            modelresol: float, optional
                An argument for model based map calculation using REFMAC.
                Resolution to calculate model based map. If not specified, an FSC
                based cutoff will be used.
            bfac: float, optional
                An overall B-factor for model map. Default is 0.0
            lig: bool, optional
                An argument for model based map calculation using REFMAC.
                Set True, if there is a ligand in the model, but no description.
                Default is True.
            lgf: string, optional
                An argument for model based map calculation using REFMAC.
                Ligand description file (cif).

        Outputs:
            fsc_list: list
                List of FSCs is returned. If len(fsc_list) is 4,
                FSC lables are as follows:
                0 - half maps FSC
                1 - half1map - model1 FSC
                2 - half2map - model1 FSC
                3 - fullmap-fullmodel FSC
                If len(fsc_list) is 2, only 0 and 3 contains.
            Outputs FSCs in allmap_fsc_modelvsmap.eps
    """
    from emda.ext import map_fsc

    fsc_list = map_fsc.map_model_fsc(
        half1_map=half1map,
        half2_map=half2map,
        modelf_pdb=modelfpdb,
        bfac=bfac,
        lig=lig,
        model1_pdb=model1pdb,
        mask_map=mask,
        model_resol=modelresol,
        lgf=lgf,
    )
    return fsc_list


def mapmodel_fsc(
    map1,
    model,
    fobj,
    bfac=0.0,
    modelresol=5.0,
    lig=True,
    phaserand=False,
    mask=None,
    lgf=None,
):
    from emda.ext import map_fsc

    res_arr, bin_fsc = map_fsc.fsc_mapmodel(
        map1=map1,
        model=model,
        model_resol=modelresol,
        bfac=bfac,
        lig=lig,
        mask_map=mask,
        lgf=lgf,
        phaserand=phaserand,
        fobj=fobj,
    )
    return res_arr, bin_fsc


def difference_map(maplist, diffmapres=3.0, ncy=100, mode="norm", fit=False, 
        usehalfmaps=False, usecom=False, fitres=None, masklist=None, 
        rot=None, axr=None, trans=None):
    """Calculates difference map.

    Arguments:
        Inputs:
            maplist: string
                List of map names to calculate difference maps.
                If combined with fit parameter, firstmap in the list
                will be taken as static/reference map. If this list
                contains coordinate file (PDB/CIF), give it in the second place.
                Always give MRC/MAP file at the beginning of the list.
                e.g:
                    [test1.mrc, test2.mrc] or
                    [test1.mrc, model1.pdb/cif]
                If combined with usehalfmaps argument, then halfmaps of the
                firstmap should be given first and those for second next.
                e.g:
                    [map1-halfmap1.mrc, map1-halfmap2.mrc, 
                     map2-halfmap1.mrc, map2-halfmap2.mrc]

            masklist: string, optional
                List of masks to apply on maps.
                All masks should be in MRC/MAP format.
                e.g:
                    [mask1.mrc, mask2.mrc]

            diffmapres: float
                Resolution to which difference map be calculated.
                If an atomic model involved, this resolution will be used
                for map calculation from coordinates

            mode: string, optional
                Different modes to scale maps. Two difference modes are supported.
                'ampli' - scale between maps is based on amplitudes .
                'norm' [Default] - normalized maps are used to calculate difference map.
                If fit is enabled, only norm mode used.

            usehalfmaps: boolean
                If employed, halfmaps are used for fitting and
                difference map calculation.
                Default is False.

            fit: boolean
                If employed, maps and superimposed before calculating
                difference map.
                Default is False.

        Outputs:
            Outputs diffence maps and initial maps after scaling in MRC format.
            Differece maps are labelled as
                emda_diffmap_m1.mrc, emda_diffmap_m2.mrc
            Scaled maps are labelled as
                emda_map1.mrc, emda_map2.mrc
    """
    from emda.ext import diffmap_with_fit

    diffmap_with_fit.difference_map(maplist=maplist, 
                                    masklist=masklist,
                                    diffmapres=diffmapres, 
                                    mode=mode, 
                                    ncy=ncy, 
                                    fit=fit, 
                                    usecom=usecom, 
                                    fitres=fitres, 
                                    usehalfmaps=usehalfmaps,
                                    rot=rot, axr=axr,
                                    trans=trans)


def applymask(mapname, maskname, outname):
    uc, arr1, origin = iotools.read_map(mapname)
    _, mask, _ = read_map(maskname)
    iotools.write_mrc(arr1 * mask, outname, uc, origin)


def scale_map2map(staticmap, map2scale, outfile):
    # this need further options
    from emda.ext.scale_maps import scale_twomaps_by_power, transfer_power

    uc, arr1, origin = iotools.read_map(staticmap)
    uc, arr2, origin = iotools.read_map(map2scale)
    f1 = np.fft.fftshift(np.fft.fftn(arr1))
    f2 = np.fft.fftshift(np.fft.fftn(arr2))
    nbin, res_arr, bin_idx = restools.get_resolution_array(uc, f1)
    scale_grid = transfer_power(
        bin_idx,
        res_arr,
        scale_twomaps_by_power(f1, f2, bin_idx=bin_idx, res_arr=res_arr),
    )
    data2write = np.real(np.fft.ifftn(np.fft.ifftshift(f2 * scale_grid)))
    iotools.write_mrc(data2write, outfile, uc, origin)
    return data2write


def bestmap(hf1name, hf2name, outfile, mode=1, knl=5, mask=None, B=None):
    from emda.ext import bestmap

    if mask is None:
        msk = 1.0
    else:
        _, msk, _ = read_map(mask)
    uc, arr1, origin = iotools.read_map(hf1name)
    uc, arr2, origin = iotools.read_map(hf2name)
    if mask:
        print("mask is not included in FSC calculation")
    f1 = np.fft.fftshift(np.fft.fftn(arr1))  # * msk))
    f2 = np.fft.fftshift(np.fft.fftn(arr2))  # * msk))
    if mode == 1:
        nbin, res_arr, bin_idx = restools.get_resolution_array(uc, f1)
        f_map = bestmap.bestmap(
            f1=f1, f2=f2, bin_idx=bin_idx, nbin=nbin, mode=mode, res_arr=res_arr, B=B)
    elif mode == 2:
        f_map = bestmap.bestmap(f1=f1, f2=f2, mode=mode, kernel_size=knl)
    data2write = np.real(np.fft.ifftn(np.fft.ifftshift(f_map))) * msk
    iotools.write_mrc(data2write, outfile, uc, origin)


def predict_fsc(hf1name, hf2name, nparticles=None, bfac=None, mask=None):
    uc, arr1, _ = iotools.read_map(hf1name)
    uc, arr2, _ = iotools.read_map(hf2name)
    if mask is not None:
        _, msk, _ = read_map(mask)
    else:
        msk = 1.0
    if nparticles is None and bfac is None:
        print("Either nparticles or bfac needs to be given!")
        exit()
    f1 = np.fft.fftshift(np.fft.fftn(arr1 * msk))
    f2 = np.fft.fftshift(np.fft.fftn(arr2 * msk))
    nbin, res_arr, bin_idx = restools.get_resolution_array(uc, f1)
    if nparticles is not None:
        bfac = None
        nparticles = 1.0 / np.asarray(nparticles, dtype="float")
        fsc_lst = fsc.predict_fsc(
            hf1=f1,
            hf2=f2,
            bin_idx=bin_idx,
            nbin=nbin,
            nparticles=nparticles,
            res_arr=res_arr,
        )
        labels = [str(i) for i in nparticles]
    if bfac is not None:
        fsc_lst = fsc.predict_fsc(
            hf1=f1, hf2=f2, bin_idx=bin_idx, nbin=nbin, bfac=bfac, res_arr=res_arr
        )
        labels = [str(i) for i in bfac]
    labels.append("reference")
    plotter.plot_nlines(
        res_arr, fsc_lst, curve_label=labels, mapname="fsc_predicted.eps"
    )
    return fsc_lst, res_arr, bin_idx, nbin


def prepare_refmac_data(
    hf1name,
    hf2name,
    outfile="fullmap.mtz",
    bfac=None,
    maskname=None,
    xmlobj=None,
    fsccutoff=None,
):
    import os
    from emda.ext import xmlclass

    xml = xmlclass.Xml()
    uc, arr1, origin = iotools.read_map(hf1name)
    uc, arr2, origin = iotools.read_map(hf2name)
    xml.map1path = os.path.abspath(hf1name)
    xml.map2path = os.path.abspath(hf2name)
    if maskname is None:
        msk = 1.0
    if maskname is not None:
        _, msk, _ = read_map(maskname)
        if arr1.shape != msk.shape:
            print("mask dim and map dim do not match!")
            print("map dim:", arr1.shape, "mask dim:", msk.shape)
            exit()
    if bfac is None:
        bfac = np.array([0.0], dtype="float")
    else:
        bfac = np.asarray(bfac, dtype="float")
        bfac = np.insert(bfac, 0, 0.0)
    hf1 = np.fft.fftshift(np.fft.fftn(arr1 * msk))
    hf2 = np.fft.fftshift(np.fft.fftn(arr2 * msk))
    # stats from half maps
    nx, ny, nz = hf1.shape
    maxbin = np.amax(np.array([nx // 2, ny // 2, nz // 2]))
    nbin, res_arr, bin_idx, sgrid = fcodes_fast.resolution_grid(
        uc, debug_mode, maxbin, nx, ny, nz
    )
    res_arr = res_arr[:nbin]
    _, _, nvar, svar, tvar, binfsc, bincount = fcodes_fast.calc_fsc_using_halfmaps(
        hf1, hf2, bin_idx, nbin, debug_mode, hf1.shape[0], hf1.shape[1], hf1.shape[2]
    )
    xml.res_arr = res_arr
    xml.fsc = binfsc
    xml.outmap = outfile
    # determine resolution
    bin_fsc = binfsc[binfsc > 0.1]
    if len(bin_fsc) > 0:
        if fsccutoff is None:
            fsccutoff = 0.5
        dist500 = np.sqrt((bin_fsc - float(fsccutoff)) ** 2)
        dist143 = np.sqrt((bin_fsc - 0.143) ** 2)
        xml.fsccutoff1 = float(0.143)
        xml.fsccutoff2 = float(fsccutoff)
        xml.mapres1 = res_arr[np.argmin(dist143)]
        xml.mapres2 = res_arr[np.argmin(dist500)]
    xml.write_xml()

    tdata = open("table_variances.txt", "w")
    tdata.write("halfmap1 file: %s\n" % os.path.abspath(hf1name))
    tdata.write("halfmap2 file: %s\n" % os.path.abspath(hf2name))
    tdata.write("\n")
    tdata.write("bin # \n")
    tdata.write("resolution (Ang.) \n")
    tdata.write("signal variance \n")
    tdata.write("noise variance \n")
    tdata.write("total variance \n")
    tdata.write("halfmap fsc \n")
    tdata.write("# reflx \n")
    for i in range(len(res_arr)):
        sv = svar[i]
        nv = nvar[i]
        tv = tvar[i]
        fsc = binfsc[i]
        nfc = bincount[i]
        tdata.write(
            "{:-3d} {:-6.2f} {:-14.4f} {:-14.4f} {:-14.4f} {:-14.4f} {:-10d}\n".format(
                i, res_arr[i], sv, nv, tv, fsc, nfc
            )
        )
    print("Bin    Resolution     FSC")
    for i in range(len(res_arr)):
        print("{:5d} {:6.2f} {:14.4f}".format(i, res_arr[i], binfsc[i]))
    # output mtz file
    iotools.write_3d2mtz_refmac(
        uc, sgrid, (hf1 + hf2) / 2.0, (hf1 - hf2), bfac, outfile=outfile
    )


def overall_cc(map1name, map2name, space="real", resol=5, maskname=None):
    """Calculates overall correlation coefficient (CC) between two density maps

    This method calculates overall CC in real or Fourier space. First two inputs can be
    both maps or one of them is map and the other is an atomic model (pdb/ent/cif).
    If one of them is an atomi model, the resolution (Angstroms) should be given.
    Default resolution is 5 A. 

    Args:
        map1name (string): Name of map1
        map2name (string): Name of map1
        space (str, optional): CC calculation in real/Fourier space. Defaults to "real".
        resol (int, optional): [description]. Defaults to 5.
        maskname ([type], optional): [description]. Defaults to None.

    Returns:
        occ: Overall CC corrected to fullmap
        hocc: Overall CC between half maps
    """
    from emda.ext import cc

    data_found = False
    if map1name.endswith((".mrc", ".map")):
        uc, arr1, origin = iotools.read_map(map1name)
        _, arr2, _ = get_data(
            struct=map2name, resol=resol, uc=uc, dim=arr1.shape, maporigin=origin
        )
    elif map1name.endswith((".pdb", ".ent", ".cif")):
        if map2name.endswith((".mrc", ".map")):
            uc, arr2, origin = iotools.read_map(map2name)
            _, arr1, _ = get_data(
                struct=map1name, resol=resol, uc=uc, dim=arr2.shape, maporigin=origin
            )
    #else:
    #    uc, arr1, origin = get_data(struct=map1name, resol=resol)
    #    _, arr1, _ = get_data(
    #        struct=map2name, resol=resol, uc=uc, dim=arr1.shape, maporigin=origin
    #    )
    if maskname is not None:
        uc, msk, origin = read_map(maskname)
        msk = msk * (msk > 0.0)
    else:
        msk = None
    if space == "fourier":
        print("Overall CC calculation in Fourier space")
        if msk is not None:
            f1 = np.fft.fftn(arr1 * msk)
            f2 = np.fft.fftn(arr2 * msk)
        else:
            f1 = np.fft.fftn(arr1)
            f2 = np.fft.fftn(arr2)
        occ, hocc = cc.cc_overall_fouriersp(f1=f1, f2=f2)
        print("Overall Correlation in Fourier space= ", occ)
    else:
        print("Overall CC calculation in Real/Image space")
        occ, hocc = cc.cc_overall_realsp(map1=arr1, map2=arr2, mask=msk)
        print("Overall Correlation in real space= ", occ)
    return occ, hocc


def mirror_map(mapname):
    """Invert map

    This method invert the map at its center-of-mass

    Args:
        mapname (string): Name of the map
    """
    # gives the inverted copy of the map
    uc, arr, origin = iotools.read_map(mapname)
    com = center_of_mass_density(arr)
    data = np.real(np.fft.ifftn(np.conjugate(np.fft.fftn(arr))))
    com2 = center_of_mass_density(data)
    data = shift_density(data, shift=np.subtract(com, com2))
    iotools.write_mrc(data, "mirror.mrc", uc, origin)


def model2map(
    modelxyz, dim, resol, cell, bfac=None, maporigin=None, ligfile=None, outputpath=None, shift_to_boxcenter=False,
):
    """Calculates EM map from atomic coordinates using REFMAC5

    Args:
        modelxyz (string): Name of the coordinate file (.cif/.pdb)
        dim (list): Map dimensions [nx, ny, nz] as a list of integers
        resol (float): Requested resolution for density calculation in Angstroms.
        cell (list): Cell parameters a, b and c as floats
        maporigin (list, optional): Location of the first column (nxstart), 
            row (nystart), section (nzstart) of the unit cell. Defaults to [0, 0, 0].
        ligfile (string, optional): Name of the ligand description file. Defaults to None.
        outputpath (string, optional): Path for auxilliary files. Defaults to current
            working directory.
        bfac(float, optional): Parameter for refmac. Set all atomic B values to bfac
            when it is positive. Default to None.
        shift_to_boxcenter (bool, optional): This parameter is useful if the calculated
            map from the model needs to be placed at the center of the box. Default to
            False. Also, the shifted model will be written to outputpath directory.

    Returns:
        float ndarray: calculated model-based density array
    """
    import gemmi as gm
    import shutil

    # print parameters
    print('Requested resolution (A): ', resol)
    print('Requested sampling: ', dim)
    print('Cell [a, b, c]: ', cell)
    if outputpath is None:
        outputpath = os.getcwd()
    outputpath = os.path.join(outputpath, 'emda_refmacfiles/')
    print('outputpath: ', outputpath)
    # make director for files for refmac run
    if os.path.exists(outputpath):
        shutil.rmtree(outputpath)
    os.mkdir(outputpath) 
    # check for valid sampling:
    for i in range(3):
        if dim[i] % 2 != 0:
            dim[i] += 1
    # check for minimum sampling
    min_pix_size = resol / 2  # in Angstrom
    min_dim = np.asarray(cell[:3], dtype="float") / min_pix_size
    min_dim = np.ceil(min_dim).astype(int)
    for i in range(3):
        if min_dim[i] % 2 != 0:
            min_dim += 1
        if min_dim[0] > dim[0]:
            print("Requested dims: ", dim)
            print("Minimum dims needed (for requested resolution): ", min_dim)
            print("!!! Please lower the requested resolution or increase the grid dimensions !!!")
            raise SystemExit()
    # replace/add cell and write model.cif
    if shift_to_boxcenter:
        from emda.core.modeltools import shift_to_origin,shift_model
        doc = shift_to_origin(modelxyz)
        doc.write_file(outputpath+"model1.cif")
        modelxyz = outputpath+"model1.cif"
    # run refmac using model.cif just created
    a, b, c = cell[:3]
    structure = gm.read_structure(modelxyz)
    structure.cell.set(a, b, c, 90.0, 90.0, 90.0)
    structure.spacegroup_hm = "P 1"
    structure.make_mmcif_document().write_file(outputpath+"model.cif")
    # run refmac using model.cif just created
    iotools.run_refmac_sfcalc(filename=outputpath+"model.cif", 
                              resol=resol, 
                              ligfile=ligfile,
                              bfac=bfac)
    modelmap = maptools.mtz2map(outputpath+"sfcalc_from_crd.mtz", dim)
    if shift_to_boxcenter:
        maporigin = None # no origin shift allowed
        modelmap = np.fft.fftshift(modelmap) #bring modelmap to boxcenter
        # shift model to boxcenter
        doc = shift_model(mmcif_file=outputpath+"model.cif", shift=[a/2, b/2, c/2])
        doc.write_file(outputpath+"emda_shifted_model.cif")
    if maporigin is None:
        maporigin = [0, 0, 0]
    else:
        shift_z = modelmap.shape[0] - abs(maporigin[2])
        shift_y = modelmap.shape[1] - abs(maporigin[1])
        shift_x = modelmap.shape[2] - abs(maporigin[0])
        modelmap = np.roll(
            np.roll(np.roll(modelmap, -shift_z, axis=0), -shift_y, axis=1),
            -shift_x,
            axis=2,
        )
    return modelmap


def model2map_gm(modelxyz, resol, dim, cell, maporigin=None, outputpath=None, shift_to_boxcenter=False,):
    import gemmi, shutil
    from servalcat.utils.model import calc_fc_fft

    if outputpath is None:
        outputpath = os.getcwd()
    outputpath = os.path.join(outputpath, 'emda_gemmifiles/')
    print('outputpath: ', outputpath)
    # make director for files for refmac run
    if os.path.exists(outputpath):
        shutil.rmtree(outputpath)
    os.mkdir(outputpath) 
    # check for valid sampling:
    for i in range(3):
        if dim[i] % 2 != 0:
            dim[i] += 1
    # check for minimum sampling
    min_pix_size = resol / 2  # in Angstrom
    min_dim = np.asarray(cell[:3], dtype="float") / min_pix_size
    min_dim = np.ceil(min_dim).astype(int)
    for i in range(3):
        if min_dim[i] % 2 != 0:
            min_dim += 1
        if min_dim[0] > dim[0]:
            print("Requested dims: ", dim)
            print("Minimum dims needed (for requested resolution): ", min_dim)
            print("!!! Please lower the requested resolution or increase the grid dimensions !!!")
            raise SystemExit()
    if shift_to_boxcenter:
        from emda.core.modeltools import shift_to_origin,shift_model
        doc = shift_to_origin(modelxyz)
        doc.write_file(outputpath+"model1.cif")
        modelxyz = outputpath+"model1.cif"
    a, b, c = cell[:3]
    st = gemmi.read_structure(modelxyz)
    st.spacegroup_hm = "P 1"
    st.cell.set(a, b, c, 90., 90., 90.)
    st.make_mmcif_document().write_file(outputpath+"model.cif")
    asu_data = calc_fc_fft(st=st, 
                           d_min=resol, 
                           source='electron', 
                           mott_bethe=True)

    """ st = gemmi.read_structure(modelxyz)
    st.spacegroup_hm = "P 1"
    st.cell.set(cell[0], cell[1], cell[2], 90., 90., 90.)
    dc = gemmi.DensityCalculatorX()
    dc.d_min = resol
    dc.rate = sample_rate = 1.5  # default
    dc.cutoff = 1.e-5  # default
    if bfac is None:
        # code inherited from servalcat. Thanks Keitaro.
        grid = resol/2/sample_rate
        b_min = min((cra.atom.b_iso for cra in st[0].all()))
        b_need = grid**2*8*numpy.pi**2/1.1 # Refmac's way
        b_add = b_need - b_min
        blur = max(0, bfac) # negative blur may cause non-positive definite in case of anisotropic Bs
        #
    dc.blur = blur
    dc.addends.subtract_z()
    dc.set_grid_cell_and_spacegroup(st)
    dc.put_model_density_on_grid(st[0])
    grid = gemmi.transform_map_to_f_phi(dc.grid)
    asu_data = grid.prepare_asu_data(
        dmin=resol, mott_bethe=True, unblur=dc.blur) """
    griddata = asu_data.get_f_phi_on_grid(
        asu_data.get_size_for_hkl(min_size=dim))
    griddata_np = (np.array(griddata, copy=False)).transpose()
    modelmap = (np.fft.ifftn(np.conjugate(griddata_np))).real
    if shift_to_boxcenter:
        maporigin = None # no origin shift allowed
        modelmap = np.fft.fftshift(modelmap) #bring modelmap to boxcenter
        # shift model to boxcenter
        doc = shift_model(mmcif_file=outputpath+"model.cif", shift=[a/2, b/2, c/2])
        doc.write_file(outputpath+"emda_shifted_model.cif")
    if maporigin is None:
        maporigin = [0, 0, 0]
    else:
        shift_z = modelmap.shape[0] - abs(maporigin[2])
        shift_y = modelmap.shape[1] - abs(maporigin[1])
        shift_x = modelmap.shape[2] - abs(maporigin[0])
        # print(shift_z, shift_y, shift_x)
        modelmap = np.roll(
            np.roll(np.roll(modelmap, -shift_z, axis=0), -shift_y, axis=1),
            -shift_x,
            axis=2,
        )
    return modelmap


def read_atomsf(atom, fpath=None):
    import subprocess

    if fpath is None:
        CMD = "echo $CLIBD"
        p = subprocess.Popen(
            CMD, stdout=subprocess.PIPE, shell=True, executable="/bin/bash"
        )
        list_of_strings = [
            x.decode("utf-8").rstrip("\n") for x in iter(p.stdout.readlines())
        ]
        fpath = list_of_strings[0] + "/atomsf.lib"
    z, ne, a, b, ier = iotools.read_atomsf(atom, fpath=fpath)
    return z, ne, a, b, ier


def compositemap(maps, masks):
    from emda.ext import composite

    composite.main(mapslist=maps, masklist=masks)


def mapmagnification(maplist, rmap, resol=4.0, masklist=None):
    """Refine magnification of a set of maps relative to a reference map

    This method refines the magnification of a set of maps relative to
    a given reference map.

    Args:
        maplist (string): List of maps whose magnifications are to be refined
        rmap (string): Name of the reference map
        masklist (string, optional): List of mask maps to apply on the maps.
        resol (float, optional): Data resolution to use in the magnification refinement.
                       This cutoff is used for all maps. Default to 4 A.
    """
    from emda.ext import magnification

    maplist.insert(0, rmap)
    magnification.main(maplist=maplist, masklist=masklist, resol=resol)


def set_dim_even(x):
    """Sets all dimentions even

    This function accepts 3D numpy array and sets its all 3 dims even

    Arguments:
        Inputs:
            x: 3D numpy array

        Outputs:
            x: 3D numpy array with all dims are even
    """
    if x.shape[0] % 2 != 0:
        xshape = list(x.shape)
        xshape[0] = xshape[0] + 1
        xshape[1] = xshape[1] + 1
        xshape[2] = xshape[2] + 1
        temp = np.zeros(xshape, x.dtype)
        temp[:-1, :-1, :-1] = x
        x = temp
    return x


def set_dim_equal(x):
    """Sets all dimentions equal and even

    This function accepts 3D numpy array and sets its all 3 dims even and equal

    Arguments:
        Inputs:
            x: 3D numpy array

        Outputs:
            x: 3D numpy array with all dims are even and equal
    """
    xshape = list(x.shape)
    maxdim = max(xshape)
    if maxdim % 2 != 0:
        maxdim = maxdim + 1
    temp = np.zeros((maxdim, maxdim, maxdim), dtype=x.dtype)
    temp[0: xshape[0], 0: xshape[1], 0: xshape[2]] = x
    x = temp
    return x


def center_of_mass_density(arr):
    """Returns the center of mass of 3D density array.

    This function accepts density as 3D numpy array and caclulates the
    center-of-mass.

    Arguments:
        Inputs:
            arr: density as 3D numpy array

        Outputs:
            com: tuple, center-of-mass (x, y, z)
    """
    from scipy import ndimage

    return ndimage.measurements.center_of_mass(arr * (arr >= 0.0))


def shift_density(arr, shift):
    """Returns a shifted copy of the input array.

    Shift the array using spline interpolation (order=3). Same as Scipy
    implementation.

    Arguments:
        Inputs:
            arr: density as 3D numpy array
            shift: sequence. The shifts along the axes.

        Outputs:
            shifted_arr: ndarray. Shifted array
    """
    from scipy import ndimage

    # return ndimage.interpolation.shift(arr, shift)
    return ndimage.shift(arr, shift, mode="wrap")


def rotate_density(arr, rotmat, interp="linear"):
    """Returns a rotated array of density

    Rotates the array of density using inperpolation.

    Arguments:
        Inputs:
            arr: density as 3D numpy array
            rotmat: 3 x 3 rotation matrix as 2D numpy array.
            interp: string.
                    Type of interpolation to use: cubic or linear.
                    Default is linear

        Outputs:
            rotated_arr: ndarray. Rotated array.
    """
    import fcodes_fast as fcodes

    nx, ny, nz = arr.shape
    if interp == "cubic":
        arr = arr.transpose()
        if arr.ndim == 3:
            arr = np.expand_dims(arr, axis=3)
            arr2 = fcodes.tricubic_map(rotmat.transpose(), arr, 1, 1, nx, ny, nz)[
                :, :, :, 0
            ]
        return arr2
    else:
        return fcodes.trilinear_map(rotmat.transpose(), arr, debug_mode, nx, ny, nz)


def get_dim(model, shiftmodel="new1.cif"):
    """Returns the box dimension to put the modelmap in.

    Determines the dimension of the box for the model based map.

    Arguments:
        Inputs:
            model:  atomic model as .pdb/.cif
            shiftmodel: name for COM shifted model, optional.
                    Default name - new1.cif.

        Outputs:
            dim: integer, dimension of the box.
    """
    import gemmi

    st = gemmi.read_structure(model)
    model = st[0]
    com = model.calculate_center_of_mass()
    print(com)
    xc = []
    yc = []
    zc = []
    for cra in model.all():
        cra.atom.pos.x -= com.x
        cra.atom.pos.y -= com.y
        cra.atom.pos.z -= com.z
        xc.append(cra.atom.pos.x)
        yc.append(cra.atom.pos.y)
        zc.append(cra.atom.pos.z)
    st.spacegroup_hm = "P 1"
    st.make_mmcif_document().write_file(shiftmodel)
    xc_np = np.asarray(xc)
    yc_np = np.asarray(yc)
    zc_np = np.asarray(zc)
    distances = np.sqrt(np.power(xc_np, 2) +
                        np.power(yc_np, 2) + np.power(zc_np, 2))
    dim1 = 2 + (int(np.max(distances)) + 1) * 2
    return dim1


def fetch_data(emdbidlist, alldata=False):
    from emda.ext import downmap

    downmap.main(emdbidlist, alldata=alldata)


""" def symaxis_refine(maplist, mapoutvar=False, emdbidlist=None, reslist=None):
    from emda.ext import refine_symaxis

    (
        emdcode_list,
        initial_axes_list,
        final_axes_list,
        fold_list,
        avgfsc_list,
    ) = refine_symaxis.main(maplist=maplist, emdbidlist=emdbidlist, mapoutvar=mapoutvar, reslist=reslist)
    return emdcode_list, initial_axes_list, final_axes_list, fold_list, avgfsc_list """


def get_map_pointgroup(maplist, reslist, use_peakheight=True, peak_cutoff=0.8,
                   use_fsc=False, fsc_cutoff=0.7, ang_tol=5.0, emdlist=None,
                   fobj=None):
    """Returns the point group of map.

    This function determines the point group of an EM map using ProSHADE and EMDA.

    Arguments:
        Inputs:
            maplist: list of strings
            List of map names in mrc/map format.

            reslist: float
            List of map resolutions

            use_peakheight: bool
            ProSHADE peaklist is used to decide the point group. Default option.

            peak_cutoff: float
            Cutoff for peak heights in the peaklist. Default is 0.8.
            However, if the highest peak is lower than this threshold then the
            cutoff is chosen such that highest - 0.1.

            use_fsc: bool
            If true, FSC is used in place of Proshade peak list to decide the point group.

            fsc_cutoff: float
            Cutoff for FSC. Default is 0.7.

            ang_tol: float
            Tolerence for angle between two axes in the proshade axes list.
            Default is 5.0 degrees.

            emdlist: list of strings
            EMDB-id list to keep the correspondence in results.


        Outputs:
            pglist: list of strings
            Point group list decided by EMDA

            ppglist: list of strings
            Point group kist decided by ProSHADE

    """
    from emda.ext.sym.symmetry import get_pointgroup

    pglist, ppglist = get_pointgroup(
        maplist=maplist,
        reslist=reslist,
        use_peakheight=use_peakheight,
        peak_cutoff=peak_cutoff,
        use_fsc=use_fsc,
        fsc_cutoff=fsc_cutoff,
        ang_tol=ang_tol,
        emdlist=emdlist,
        fobj=fobj)
    return pglist, ppglist


def symmetry_average(maplist, reslist, use_peakheight=True, peak_cutoff=0.8,
                   use_fsc=False, fsc_cutoff=0.7, ang_tol=5.0, pglist=None,
                   emdlist=None, fobj=None):
    """Returns symmetry averaged map.

    This function does three difference things:
        1. Determines the point group of a map using ProSHADE.
        2. Identify group generators and refine them.
        3. Generate the full finite point group operators and use them to average the map.
    If point groups are supplied, then they are used for symmetry averaging in
    C, D and O groups. Averaging in T and I groups is done using operators
    from refined axes.
    This function can be used to average maps, whose point groups
    and symmetry operators are not known a priori.
    NOTE: Current axis refinement for T and I groups may include numerical
    inaccuracies and so the symmetry averaged map may not be the optimal.

    Arguments:
        Inputs:
            maplist: list of strings
            List of map names in mrc/map format.

            reslist: float
            List of map resolutions

            use_peakheight: bool
            ProSHADE peaklist is used to decide the point group. Default option.

            peak_cutoff: float
            Cutoff for peak heights in the peaklist. Default is 0.8.
            However, if the highest peak is lower than this threshold then the
            cutoff is chosen such that highest - 0.1.

            use_fsc: bool
            If true, FSC is used in place of Proshade peak list to decide the point group.

            fsc_cutoff: float
            Cutoff for FSC. Default is 0.7.

            ang_tol: float
            Tolerence for angle between two axes in the proshade axes list.
            Default is 5.0 degrees.

            emdlist: list of strings
            EMDB-id list to keep the correspondence in results.

            pglist: list of strings
            List of point groups

        Outputs:
            symavgmaplist: list of maps
            List of symmetry averaged maps
    """
    from emda.ext.sym.symmetry import symmetrise_map

    symavgmaplist = symmetrise_map(
        maplist=maplist,
        reslist=reslist,
        use_peakheight=use_peakheight,
        peak_cutoff=peak_cutoff,
        use_fsc=use_fsc,
        fsc_cutoff=fsc_cutoff,
        ang_tol=ang_tol,
        pglist=pglist,
        emdlist=emdlist)
    return symavgmaplist


def symmetry_average_using_ops(imap, ops, outmapname=None):
    """Returns symmetry averaged map using given operators.

    This function can be used to average a map by a given set of symmetry operators.

    Arguments:
        Inputs:
            imap:  3D-EM map in mrc/map format.
            ops: List of symmetry operators to be applied on the map.
                 List should not contain the identity.
                 An operator must be a 3x3 matrix.

        Outputs:
            Results in a list: [symmetry-averaged-density, unit-cell, map-origin]
    """
    from emda.ext.sym.symmetrize_map import symmetrize_map_using_ops

    return symmetrize_map_using_ops(imap=imap, ops=ops, outmapname=outmapname)


def rebox_mapmodel(maplist, masklist, modellist=None):
    """Reboxes a map by a mask.

    This function reboxes a map by a mask into a smaller cube. Size of the box 
    is determine by the size of te mask + 10 voxels in each dimetion. 
    It can accept list arguments if there's more than one map and mask. 
    Also, a model can be given and then that will also be reboxed.

    Args:
        maplist (string): List of map names to be reboxed
        masklist (string): List of mask names to be used in reboxing
        modellist (string, optional): List of model names (pdb/cif) to be reboxed. 
                        Defaults to None.
        
    Output:
        It outputs reboxed maps, models (if supplied).
    """
    from emda.ext.rebox_map import rebox_maps_and_models

    rebox_maps_and_models(maplist=maplist, masklist=masklist, modellist=modellist)