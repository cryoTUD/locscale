# emda test exhaust
# each method in emda_methods is tested against test data

from __future__ import absolute_import, division, print_function, unicode_literals
import fcodes_fast
from emda import emda_methods as em

map1name = (
    "/Users/ranganaw/MRC/REFMAC/haemoglobin/EMD-3651/other/emd_3651_half_map_1.map"
)
map2name = (
    "/Users/ranganaw/MRC/REFMAC/haemoglobin/EMD-3651/other/emd_3651_half_map_2.map"
)
modelf = "/Users/ranganaw/MRC/REFMAC/haemoglobin/EMD-3651/refined.pdb"

# Tests
def main():
    uc, arr1, orig = em.read_map(map1name)
    uc, arr2, orig = em.read_map(map2name)
    print("reading test - PASS")

    em.write_mrc(arr1, "test.mrc", uc, map_origin=orig)
    print("writing test - PASS")

    new_arr = em.resample_data(
        curnt_pix=[1.0,1.0, 1.0], targt_pix=[1.0,1.0, 1.0], targt_dim=[100, 100, 100], arr=arr1
    )
    print("resampling test - PASS")

    resol = em.estimate_map_resol(hfmap1name=map1name, hfmap2name=map2name)
    print("estimate resolution test - PASS")

    res_arr, power_spectrum = em.get_map_power(map1name)
    print("calculate power spectrum test - PASS")

    all_mapout = em.apply_bfactor_to_map(map1name, bf_arr=[0.0], mapout=True)
    print("applying B factor to map test - PASS")

    em.map2mtz(map1name, mtzname="test.mtz")
    print("writing map to mtz test - PASS")

    mtz2map = em.mtz2map(mtzname="test.mtz", map_size=arr1.shape)
    print("writing mtz to map test - PASS")

    fmap1, map1 = em.lowpass_map(uc, arr1, resol, filter="ideal", order=4)
    print("lowpass filtering test - PASS")

    fullmap = em.half2full(
        half1name=map1name, half2name=map2name, outfile="fullmap.mrc"
    )
    print("fullmap from halfmaps test - PASS")

    transformedmap = em.map_transform(
        mapname=map1name, tra=[0, 0, 0], rot=0, axr=[1, 0, 0], outname="transformedmap.mrc"
    )
    print("apply transformation on map test - PASS")

    res_arr, bin_fsc = em.halfmap_fsc(half1name=map1name, half2name=map2name)
    print("FSC from halfmaps test - PASS")

    res_arr, noisevar, signalvar = em.get_variance(
        half1name=map1name, half2name=map2name
    )
    print("calculating variances from halfdata test - PASS")

    res_arr, bin_fsc = em.twomap_fsc(map1name=map1name, map2name=map2name)
    print("FSC between two maps test - PASS")

    mask = em.mask_from_halfmaps(uc=uc, half1=arr1, half2=arr2, radius=5)
    em.write_mrc(mask, "mask.mrc", uc, orig)
    print("ccmask test - PASS")

    mask = em.mask_from_map(arr=arr1, kern=5, uc=uc)
    print("mapmask test - PASS")

    em.overlay_maps(maplist=[map1name, map2name], masklist=["mask.mrc", "mask.mrc"])
    print("map overlay test - PASS")

    em.realsp_correlation(half1map=map1name, half2map=map2name, kernel_size=5)
    print("local correlation rcc test - PASS")

    em.fouriersp_correlation(half1_map=map1name, half2_map=map2name, kernel_size=5)
    print("local correlation fcc test - PASS")

    modelmap = em.model2map(modelxyz=modelf, dim=arr1.shape, resol=resol, cell=uc)
    em.write_mrc(modelmap, "modelmap.mrc", uc, orig)
    print("model2map test - PASS")

    em.map_model_validate(
        half1map=map1name,
        half2map=map2name,
        modelfpdb=modelf,
        model1pdb=modelf
    )

    em.difference_map(
        maplist=[map1name, "modelmap.mrc"],
        masklist=["mask.mrc", "mask.mrc"],
        diffmapres=resol,
        mode="ampli"
    )


if __name__ == "__main__":
    main()
