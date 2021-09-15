"""
Author: "Rangana Warshamanage, Garib N. Murshudov"
MRC Laboratory of Molecular Biology
    
This software is released under the
Mozilla Public License, version 2.0; see LICENSE.
"""

# local correlation in real space

from __future__ import absolute_import, division, print_function, unicode_literals
import numpy as np
from emda import core
import emda.emda_methods as em

class RealspaceLocalCC:
    def __init__(self):
        self.hfmap1name = None
        self.hfmap2name = None
        self.maskname = None
        self.uc = None
        self.origin = None
        self.kern_rad = 5
        self.norm=False
        self.model=None
        self.model_resol=None
        self.lgf=None
        self.arr1 = None
        self.arr2 = None
        self.mask = None
        self.hfmapcc = None
        self.fullmapcc = None
        self.mapmodelcc = None
        self.truemapmodelcc = None

    def check_inputs(self):
        if self.model is not None:
            if self.model.endswith((".pdb", ".ent", ".cif")):
                if self.model_resol is None:
                    raise SystemExit(
                        "Please input resolution for model-based map!")

    def get_fullmapcor(self):
        import numpy.ma as ma

        hfmap_corr_ma = ma.masked_less_equal(self.hfmapcc, 0.0)
        ccf_ma = 2 * hfmap_corr_ma / (1.0 + hfmap_corr_ma)
        self.mapmodelcc = ccf_ma.filled(0.0)

    def truemap_model_cc(self):
        import numpy.ma as ma

        # To prevent large numbers in truemapmodelcc
        mapmodelcc_ma = ma.masked_less_equal(self.mapmodelcc, 0.3)
        fullmapcc_ma = ma.masked_less_equal(self.fullmapcc, 0.3)
        truemapmodelcc_ma = mapmodelcc_ma / np.sqrt(fullmapcc_ma)
        truemapmodelcc = core.iotools.mask_by_value_greater(
            truemapmodelcc_ma.filled(0.0), masking_value=1.0
        )
        self.truemapmodelcc = truemapmodelcc #* (truemapmodelcc > 0.0)

    def outputmaps(self, data, outname):
        core.iotools.write_mrc(
            mapdata=data,
            filename=outname + str(self.kern_rad) + ".mrc",
            unit_cell=self.uc,
            map_origin=self.origin,
            label=True,
        )

    def rcc(self):
        self.check_inputs()
        model = self.model
        print(
            "Calculating 3D correlation using half maps. \
             Please wait..."
        )
        uc, arr1, origin = core.iotools.read_map(self.hfmap1name)
        uc, arr2, origin = core.iotools.read_map(self.hfmap2name)
        self.uc, self.origin = uc, origin
        if self.maskname is not None:
            _, cc_mask, _ = core.iotools.read_map(self.maskname)
            cc_mask_binary = cc_mask > 0.5 # binary mask
        else:
            cc_mask = 1
            cc_mask_binary = 1
        self.mask = cc_mask
        nx, ny, nz = arr1.shape
        kern_sphere_soft = core.restools.create_soft_edged_kernel_pxl(self.kern_rad)
        f_hf1, f_hf2 = np.fft.fftn(arr1*self.mask), np.fft.fftn(arr2*self.mask)
        print("Calculating 3D correlation...")
        if self.norm:
            hf1, hf2 = np.fft.fftshift(f_hf1), np.fft.fftshift(f_hf2)
            print("Normalising maps...")
            nm = NormalizedMaps(hf1=hf1, hf2=hf2, cell=uc)
            nm.get_normdata()
            nbin, bin_idx = nm.nbin, nm.bin_idx
            normfull = nm.normfull
            em.write_mrc(nm.normmap1, "bin_normalized_halfmap1.mrc", uc, origin)
            em.write_mrc(nm.normmap2, "bin_normalized_halfmap2.mrc", uc, origin)
            em.write_mrc(normfull, "bin_normalized_fullmap.mrc", uc, origin)
            # Real space correlation maps
            print("Calculating 3D correlation using normalized maps...")
            halfmapscc = get_3d_realspcorrelation(
                half1=nm.normmap1 * cc_mask,
                half2=nm.normmap2 * cc_mask,
                kern=kern_sphere_soft,
            )
        if not self.norm:
            halfmapscc = get_3d_realspcorrelation(
                half1=arr1 * self.mask, 
                half2=arr2 * self.mask, 
                kern=kern_sphere_soft
            )
        fullmapcc = 2 * halfmapscc / (1.0 + halfmapscc)
        self.hfmapcc = halfmapscc 
        self.fullmapcc = fullmapcc
        print("Writing out correlation maps...")
        self.outputmaps(self.hfmapcc * cc_mask_binary, "rcc_halfmap_smax")
        self.outputmaps(self.fullmapcc * cc_mask_binary, "rcc_fullmap_smax")
        self.outputmaps(np.sqrt(
                        np.where(self.fullmapcc <= 0.0, 0.0, self.fullmapcc))
                        * cc_mask_binary, "rcc_fullmap_star_smax")
        # Map-model correlation
        if model is not None:
            print("\nMap-model correlation!\n")
            dim = [nx, ny, nz]
            if model.endswith((".pdb", ".ent", ".cif")):
                print("Map will be calculated up to " 
                      + str(self.model_resol) + "A")
                print("Calculating map from model using REFMAC. Please wait...")
                model_arr = calculate_modelmap(
                    modelxyz=model, 
                    dim=dim, 
                    resol=self.model_resol, 
                    uc=uc, 
                    maporigin=origin,
                    lgf=self.lgf)
            elif model.lower().endswith((".mrcs", ".mrc", ".map")):
                _, model_arr, _ = core.iotools.read_map(model)
            elif model.lower().endswith(".mtz"):
                model_arr = core.maptools.mtz2map(model, (nx, ny, nz))
            em.write_mrc(model_arr, "modelmap.mrc", uc, origin)
            if self.norm:
                print("Calculating model-map correlation...\n")
                normmodel = normalized_modelmap(
                                modelmap=model_arr, 
                                fsc_grid_str=nm.fsc_grid_str, 
                                bin_idx=bin_idx, 
                                nbin=nbin
                            )
                em.write_mrc(normmodel, "bin_normalized_modelmap.mrc", uc, origin)
                print("Calculating model-map correlation...\n")
                mapmodelcc = get_3d_realspcorrelation(
                    half1=normfull * cc_mask,
                    half2=normmodel * cc_mask,
                    kern=kern_sphere_soft,
                )
            if not self.norm:
                f_fullmap = (f_hf1 + f_hf2) / 2.0
                fullmap = np.real(np.fft.ifftn(f_fullmap))
                print("Calculating model-map correlation...\n")
                mapmodelcc= get_3d_realspcorrelation(fullmap, model_arr, kern_sphere_soft)
            self.mapmodelcc = mapmodelcc
            self.outputmaps(self.mapmodelcc * cc_mask_binary, "rcc_mapmodel_smax")
            print("Calculating truemap-model correlation...")
            #self.truemap_model_cc()
            #self.outputmaps(self.truemapmodelcc * cc_mask_binary, "rcc_truemapmodel_smax")
            print("Map-model correlation calculated! Maps were writted!")


def get_3d_realspcorrelation(half1, half2, kern, mask=None):
    import scipy.signal
    from scipy.stats import mode

    loc3_A = scipy.signal.fftconvolve(half1, kern, "same")
    loc3_A2 = scipy.signal.fftconvolve(half1 * half1, kern, "same")
    loc3_B = scipy.signal.fftconvolve(half2, kern, "same")
    loc3_B2 = scipy.signal.fftconvolve(half2 * half2, kern, "same")
    loc3_AB = scipy.signal.fftconvolve(half1 * half2, kern, "same")
    cov3_AB = loc3_AB - loc3_A * loc3_B
    var3_A = loc3_A2 - loc3_A ** 2
    var3_B = loc3_B2 - loc3_B ** 2
    # regularization
    #reg_a = mode(var3_A)[0][0][0][0] / 100
    #reg_b = mode(var3_B)[0][0][0][0] / 100
    reg_a = np.max(var3_A) / 1000
    reg_b = np.max(var3_B) / 1000
    var3_A = np.where(var3_A < reg_a, reg_a, var3_A)
    var3_B = np.where(var3_B < reg_b, reg_b, var3_B)
    halfmaps_cc = cov3_AB / np.sqrt(var3_A * var3_B) 
    return halfmaps_cc


def calculate_modelmap(modelxyz, dim, resol, uc, lgf=None, maporigin=None):
    modelmap = em.model2map(
        modelxyz=modelxyz,
        dim=dim,
        resol=resol,
        cell=uc,
        ligfile=lgf,
        maporigin=maporigin,
    )
    return modelmap


def make_c_B_table(s_max, dat_out=None):
    import scipy.special
    t = lambda B: np.sqrt(np.abs(B)/2)*s_max
    w = lambda z: np.exp(-z**2)*(1-1j*scipy.special.erfi(-z))
    fp = lambda x: 3/2*(np.sqrt(np.pi)/2*scipy.special.erf(t(x))-t(x)*np.exp(-t(x)**2))/t(x)**3
    fn = lambda x: 3/4*np.exp(t(x)**2)*(2*t(x)-np.sqrt(np.pi)*np.imag(w(t(x))))/t(x)**3

    Bn = np.arange(-300,0)
    cn = fn(Bn)
    Bp = np.arange(1, 600)
    cp = fp(Bp)

    Bs = np.concatenate((Bn, [0.], Bp))
    cs = np.concatenate((cn, [1.], cp))

    order = np.argsort(cs)
    Bs, cs = Bs[order], cs[order]

    if dat_out:
      with open(dat_out, "w") as ofs:
        ofs.write("# s_max= {}\n".format(s_max))
        ofs.write("c B\n")
        for i in range(len(Bs)):
          ofs.write("{:.5e} {:.1f}\n".format(cs[i], Bs[i]))

    return Bs, cs

def bfromcc(
    half1_map,
    half2_map,
    kernel_size,
    resol,
    mask_map=None,
):
    import emda.emda_methods as em

    hf1, hf2 = None, None
    bin_idx = None
    print(
        "Calculating 3D correlation between half maps. \
            Please wait..."
    )
    uc, arr1, origin = core.iotools.read_map(half1_map)
    uc, arr2, origin = core.iotools.read_map(half2_map)

    # mask taking into account
    if mask_map is not None:
        _, cc_mask, _ = core.iotools.read_map(mask_map)
    else:
        # creating ccmask from half data
        print("Mask is not given. EMDA will generate cc_mask.mrc and be used.")
        print("Please take a look at this automatic mask. It may be suboptimal.")
        cc_mask = em.mask_from_halfmaps(uc, arr1, arr2, radius=7, thresh=0.65)
        em.write_mrc(cc_mask, "cc_mask.mrc", uc, origin)

    cc_mask = cc_mask * (cc_mask > 0.0)
    # Creating soft-edged mask
    kern_sphere_soft = core.restools.create_soft_edged_kernel_pxl(kernel_size)

    # Lowpass
    f_hf1, arr1 = em.lowpass_map(uc=uc, arr1=arr1, resol=resol)
    f_hf2, arr2 = em.lowpass_map(uc=uc, arr1=arr2, resol=resol)

    # Real space correlation maps
    print("Calculating 3D correlation...")
    halfmapscc, fullmapcc = get_3d_realspcorrelation(
        half1=arr1 * cc_mask, half2=arr2 * cc_mask, kern=kern_sphere_soft
    )

    overall_cc = np.corrcoef((arr1*cc_mask).flatten(), (arr2*cc_mask).flatten())[1,0]
    print("Overall CC=", overall_cc)

    c = (1.-overall_cc)/overall_cc * halfmapscc/(1.-halfmapscc)
    print("c range: {:.4e} to {:.4e}".format(c.min(), c.max()))
    core.iotools.write_mrc(c,'c_for_localB'+str(kernel_size)+'.mrc',uc,origin)

    table_B, table_c = make_c_B_table(1./resol, dat_out="c_B.dat")
    print("c to B table made for c= {:.4e} to {:.4e} and B= {:.1f} to {:.1f}".format(min(table_c), max(table_c), min(table_B), max(table_B)))
    B = np.interp(c, table_c, table_B)
    print("Writing out local B map")
    core.iotools.write_mrc(mapdata=B,filename='localB'+str(kernel_size)+'.mrc',unit_cell=uc,map_origin=origin)

    print("Writing out correlation maps")
    core.iotools.write_mrc(
        mapdata=halfmapscc * cc_mask,
        filename="rcc_halfmap_smax" + str(kernel_size) + ".mrc",
        unit_cell=uc,
        map_origin=origin,
        label=True,
    )


def scale_model2map(fullmap, model_arr, uc):
    from emda.ext.scale_maps import scale_twomaps_by_power, transfer_power

    f1 = np.fft.fftshift(np.fft.fftn(np.fft.fftshift(fullmap)))
    f2 = np.fft.fftshift(np.fft.fftn(np.fft.fftshift(model_arr)))
    nbin, res_arr, bin_idx = core.restools.get_resolution_array(uc, f1)
    scale_grid = transfer_power(
        bin_idx,
        res_arr,
        scale_twomaps_by_power(f1, f2, bin_idx=bin_idx, res_arr=res_arr),
    )
    scaled_model = np.real(np.fft.ifftn(np.fft.ifftshift(f2 * scale_grid)))
    return scaled_model


def mapmodel_rcc(
    fullmap,
    model,
    resol,
    kernel_size=5,
    norm=False,
    mask_map=None,
    lgf=None,
):
    try:
        print(
            "Calculating 3D correlation between map and model. \
            Please wait..."
        )
        uc, arr1, origin = core.iotools.read_map(fullmap)
        nx, ny, nz = arr1.shape
        # lowpass filter the map to given resolution
        _, arr1 = em.lowpass_map(
            uc=uc, 
            arr1=arr1, 
            resol=resol, 
            filter='butterworth')
        if mask_map is not None:
            print("Correlation is calculated using a mask.")
            print("Input map is masked before calculating correlation.")
            _, mask, _ = core.iotools.read_map(mask_map)
            mask_binary = mask > 0.
        else:
            mask = 1
            mask_binary = 1
        print("\nMap-model correlation!\n")
        dim = [nx, ny, nz]
        if model.endswith((".pdb", ".ent", ".cif")):
            print("map will be calculated up to " + str(resol) + "A")
            print("Calculating map from model using REFMAC!")
            model_arr = calculate_modelmap(
                modelxyz=model, 
                dim=dim, 
                resol=resol, 
                uc=uc, 
                lgf=lgf, 
                maporigin=origin,
            )
        elif model.lower().endswith((".mrc", ".map")):
            _, model_arr, _ = core.iotools.read_map(model)
            _, model_arr = em.lowpass_map(
                uc=uc, 
                arr1=model_arr, 
                resol=resol, 
                filter='butterworth')
        elif model.lower().endswith(".mtz"):
            model_arr = core.maptools.mtz2map(model, (nx, ny, nz))
            _, model_arr = em.lowpass_map(
                uc=uc, 
                arr1=model_arr, 
                resol=resol, 
                filter='butterworth')
        core.iotools.write_mrc(
            mapdata=model_arr, filename="modelmap.mrc", unit_cell=uc, map_origin=origin
        )
        print("Calculating model-map correlation...\n")
        kern_sphere_soft = core.restools.create_soft_edged_kernel_pxl(kernel_size)
        if norm:
            print("Normalisation needs half maps. If you have half maps, then \
                use rcc instead of mmcc. For now, local cc map-model will be \
                    calculated using un-normalized map and model.")
        mapmodelcc = get_3d_realspcorrelation(
                arr1 * mask, model_arr, kern_sphere_soft
            )
        core.iotools.write_mrc(
            mapdata=mapmodelcc * mask_binary,
            filename="mmcc_mapmodel.mrc",
            unit_cell=uc,
            map_origin=origin,
            label=True,
        )
        print("Map-model correlation calculated! Maps were written!")
    except:
        print('Input error!')


def normalized(map, bin_idx=None, nbin=None, uc=None):
    import fcodes_fast as fc

    # normalise in resol bins
    if bin_idx is None:
        if uc is not None:
            nbin, res_arr, bin_idx = core.restools.get_resolution_array(uc, map)
    hf1 = np.fft.fftshift(np.fft.fftn(map))
    binfsc, _, _, totalvar, _, eo = core.fsc.halfmaps_fsc_variance(
        hf1, hf1, bin_idx, nbin
    )
    binfsc = binfsc * (totalvar > 1e-5)
    nx, ny, nz = eo.shape
    fsc_frid = fc.read_into_grid(bin_idx, binfsc, nbin, nx, ny, nz)
    norm_map = np.real(np.fft.ifftn(np.fft.ifftshift(eo * fsc_frid)))
    return norm_map


""" def normalized_modelmap(modelmap, f_fullmap, bin_idx=None, nbin=None, uc=None):
    import fcodes_fast as fc

    # normalise in resol bins
    if bin_idx is None:
        if uc is not None:
            nbin, res_arr, bin_idx = core.restools.get_resolution_array(uc, modelmap)
    f1 = np.fft.fftshift(np.fft.fftn(modelmap))
    f_fullmap = np.fft.fftshift(f_fullmap)
    #_, _, _, _, _, eo = core.fsc.halfmaps_fsc_variance(f1, f1, bin_idx, nbin)
    #binfsc, _, _, _, _, _ = core.fsc.halfmaps_fsc_variance(f1, f_fullmap, bin_idx, nbin)
    nx, ny, nz = f1.shape
    eo, _, _, totalvar, binfsc, _ = fc.get_normalized_sf(
        f1, f_fullmap, bin_idx, nbin, 0, nx, ny, nz
    )
    e1 = eo[:, :, :, 0]
    #fsc_frid = fc.read_into_grid(bin_idx, binfsc, nbin, nx, ny, nz)
    #norm_map = np.real(np.fft.ifftn(np.fft.ifftshift(e1 * fsc_frid)))
    norm_map = np.real(np.fft.ifftn(np.fft.ifftshift(e1))) # not mulpiplying by fsc becoz no noise in the model
    return norm_map """

def normalized_modelmap(modelmap, fsc_grid_str, bin_idx=None, nbin=None, uc=None):
    import fcodes_fast as fc

    # normalise in resol bins
    if bin_idx is None:
        if uc is not None:
            nbin, res_arr, bin_idx = core.restools.get_resolution_array(uc, modelmap)
    f1 = np.fft.fftshift(np.fft.fftn(modelmap))
    nx, ny, nz = f1.shape
    eo = fc.get_normalized_sf(f1, f1, bin_idx, nbin, 0, nx, ny, nz)[0]
    norm_map = np.real(np.fft.ifftn(np.fft.ifftshift(fsc_grid_str * eo[:, :, :, 0])))
    return norm_map


def hfdata_normalized(hf1, hf2, bin_idx=None, nbin=None, uc=None):
    # normalise in resol bins
    import fcodes_fast as fc

    print("Normalising maps...")
    nm = NormalizedMaps(hf1=hf1, hf2=hf2, cell=uc)
    nm.get_normdata()
    return nm.normmap1, nm.normmap2


class NormalizedMaps:
    def __init__(self, hf1=None, hf2=None, bin_idx=None, nbin=None, cell=None):
        self.hf1 = hf1
        self.hf2 = hf2
        self.e0 = None
        self.e1 = None
        self.e2 = None
        self.normmap1 = None
        self.normmap2 = None
        self.normfull = None
        self.bin_idx = bin_idx
        self.cell = cell
        self.nbin = nbin
        self.binfsc = None
        self.res_arr = None
        self.totalvar = None
        self.fsc_grid_str = None

    def get_parameters(self):
        if self.bin_idx is None:
            if self.cell is not None:
                (
                    self.nbin,
                    self.res_arr,
                    self.bin_idx,
                ) = core.restools.get_resolution_array(self.cell, self.hf1)

    def get_normdata(self):
        import fcodes_fast as fc
        from emda.ext.mapfit.mapaverage import set_array

        if self.bin_idx is None:
            self.get_parameters()
        nx, ny, nz = self.hf1.shape
        eo, _, _, self.totalvar, self.binfsc, _ = fc.get_normalized_sf(
            self.hf1, self.hf2, self.bin_idx, self.nbin, 0, nx, ny, nz
        )
        self.e1 = eo[:, :, :, 0]
        self.e2 = eo[:, :, :, 1]
        self.e0 = eo[:, :, :, 2]
        fsc_half = set_array(arr=self.binfsc, thresh=0.01)
        fsc_ful = 2 * fsc_half / (1 + fsc_half)
        fsc_full_grid = fc.read_into_grid(self.bin_idx, fsc_ful, self.nbin, nx, ny, nz)
        fsc_grid_str = np.sqrt(fsc_full_grid)
        self.fsc_grid_str = fsc_grid_str
        self.normmap1 = np.real(np.fft.ifftn(np.fft.ifftshift(self.e1 * fsc_grid_str)))
        self.normmap2 = np.real(np.fft.ifftn(np.fft.ifftshift(self.e2 * fsc_grid_str)))
        self.normfull = np.real(np.fft.ifftn(np.fft.ifftshift(self.e0 * fsc_grid_str)))


if __name__=="__main__":
    #path = "/Users/ranganaw/MRC/REFMAC/EMD-6952/emda_test/test_mmcc/new_mmcc/"
    #fullmap = "emd_6952.map"
    #model = "modelmap.mrc"
    resol = 4.
    path = "/Users/ranganaw/MRC/REFMAC/COVID19/EMD-11203/test_rcc/"
    halfmap1 = path + "emd_11203_half1.map"
    halfmap2 = path + "emd_11203_half2.map"

    #mapmodel_rcc(fullmap, model, resol)
    mapmodel_rcc(halfmap1, halfmap2, resol)

    # calculating RCC
    #rcc = RealspaceLocalCC()
    #rcc.hfmap1name = halfmap1
    #rcc.hfmap2name = halfmap2
    #rcc.rcc()

