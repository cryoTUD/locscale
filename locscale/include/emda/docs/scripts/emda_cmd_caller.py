"""
Author: "Rangana Warshamanage, Garib N. Murshudov"
MRC Laboratory of Molecular Biology

This software is released under the
Mozilla Public License, version 2.0; see LICENSE.
"""
from __future__ import absolute_import, division, print_function, unicode_literals
import numpy as np
import argparse
import sys
import datetime
#import emda.config
#from emda.core import iotools, maptools, restools, plotter, fsc, quaternions

cmdl_parser = argparse.ArgumentParser(
    prog="emda",
    usage="%(prog)s [command] [arguments]",
    formatter_class=argparse.RawDescriptionHelpFormatter,
)

#cmdl_parser.add_argument(
#    "--version", action="version", version="%(prog)s-" + emda.config.__version__
#)

subparsers = cmdl_parser.add_subparsers(dest="command")

mapinfo = subparsers.add_parser(
    "info", description="Output basic information about the map."
)
mapinfo.add_argument("--map", required=True, help="input map")

halffsc = subparsers.add_parser(
    "halffsc", description="Calculates FSC between half-maps."
)
halffsc.add_argument("--h1", required=True, help="input map1 map")
halffsc.add_argument("--h2", required=True, help="input map2 map")
halffsc.add_argument(
    "--msk", required=False, default=None, type=str, help="input mask (mrc/map)"
)
halffsc.add_argument(
    "--out", required=False, default="table_variances.txt", help="output data table"
)
halffsc.add_argument(
    "--phaserand", action="store_true", help="use if phase randomized FSC is calculated"
)

anyfsc = subparsers.add_parser(
    "fsc", description="Calculates FSC between two maps.")
anyfsc.add_argument("--map1", required=True, help="input map1 map")
anyfsc.add_argument("--map2", required=True, help="input map2 map")

singlemapfsc = subparsers.add_parser(
    "singlemapfsc", description="Calculates FSC using neighbour average."
)
singlemapfsc.add_argument("--h1", required=True, help="input map1 map")
singlemapfsc.add_argument(
    "--knl", required=False, type=int, default=3, help="kernel radius in voxels"
)

ccmask = subparsers.add_parser(
    "ccmask", description="Generates a mask based on halfmaps correlation."
)
ccmask.add_argument("--h1", required=True, help="input halfmap1 map")
ccmask.add_argument("--h2", required=True, help="input halfmap2 map")
ccmask.add_argument(
    "--knl", required=False, type=int, default=4, help="kernel radius in voxels"
)
ccmask.add_argument(
    "--itr", required=False, type=int, default=1, help="number of dilation cycles"
)
ccmask.add_argument(
    "--dthreshold", required=False, default=None, type=float, help="threshold for density"
)

map_mask = subparsers.add_parser(
    "mapmask", description="Generate a mask from a map.")
map_mask.add_argument("--map", required=True, help="input map")
map_mask.add_argument(
    "--knl", required=False, type=int, default=5, help="kernel radius in voxels"
)
map_mask.add_argument(
    "--prb",
    required=False,
    type=float,
    default=0.99,
    help="density cutoff probability in cumulative density function",
)
map_mask.add_argument(
    "--itr", required=False, type=int, default=3, help="number of dilate iterations"
)
map_mask.add_argument(
    "--res",
    required=False,
    type=float,
    default=15.0,
    help="lowpass resolution in Angstroms",
)
map_mask.add_argument(
    "--fil",
    required=False,
    type=str,
    default="butterworth",
    help="filter type to use: ideal or butterworth",
)

model_mask = subparsers.add_parser(
    "modelmask", description="Generate a mask from an atomic model.")
model_mask.add_argument("--map", required=True, help="input map MRC/MAP")
model_mask.add_argument("--mdl", required=True,
                        help="input atomic model PDB/CIF")
model_mask.add_argument("--atmrad", required=False, default=3.0,
                        type=float, help="radius of the atomic sphere in Angstroms")

lowpass = subparsers.add_parser(
    "lowpass", description="Lowpass filter to specified resolution."
)
lowpass.add_argument("--map", required=True, help="input map (mrc/map)")
lowpass.add_argument(
    "--res", required=True, type=float, help="lowpass resolution in Angstrom"
)
lowpass.add_argument(
    "--fil",
    required=False,
    type=str,
    default="ideal",
    help="filter type to use: ideal or butterworth",
)

power = subparsers.add_parser(
    "power", description="Calculates power spectrum.")
power.add_argument("--map", required=True, help="input map (mrc/map)")

applybfac = subparsers.add_parser(
    "bfac", description="Apply a B-factor on the map.")
applybfac.add_argument("--map", required=True, help="input map (mrc/map)")
applybfac.add_argument(
    "--bfc", required=True, nargs="+", type=float, help="bfactor(s) to apply"
)
applybfac.add_argument("--out", action="store_true",
                       help="if use, writes out map")

map_resol = subparsers.add_parser(
    "resol", description="Estimates map resolution based on FSC."
)
map_resol.add_argument("--h1", required=True, help="input halfmap1 map")
map_resol.add_argument("--h2", required=True, help="input halfmap2 map")

half2full = subparsers.add_parser(
    "half2full", description="Combine two halfmaps to make the fullmap."
)
half2full.add_argument("--h1", required=True, help="input halfmap1 map")
half2full.add_argument("--h2", required=True, help="input halfmap2 map")
half2full.add_argument(
    "--out", required=False, default="fullmap.mrc", help="output map (mrc/map)"
)

conv_map2mtz = subparsers.add_parser(
    "map2mtz", description="Convert MRC/MAP to MTZ.")
conv_map2mtz.add_argument("--map", required=True, help="input map (mrc/map)")
conv_map2mtz.add_argument("--res", required=False, type=float,
                          help="resolution cutoff (A). default Nyquist")
conv_map2mtz.add_argument(
    "--out", required=False, default="map2mtz.mtz", help="output map (mtz)"
)

conv_map2mtzful = subparsers.add_parser(
    "map2mtzfull", description="Convert MRC/MAP to MTZ using half maps."
)
conv_map2mtzful.add_argument(
    "--h1", required=True, help="input hfmap1 (mrc/map)")
conv_map2mtzful.add_argument(
    "--h2", required=True, help="input hfmap2 (mrc/map)")
conv_map2mtzful.add_argument(
    "--out", required=False, default="map2mtzfull.mtz", help="output map (mtz)"
)

transform_map = subparsers.add_parser(
    "transform", description="Apply a transformation on the map."
)
transform_map.add_argument("--map", required=True, help="input map (mrc/map)")
transform_map.add_argument(
    "--tra",
    required=False,
    default=[0.0, 0.0, 0.0],
    nargs="+",
    type=float,
    help="translation vec. in Angstrom. eg 1.0 0.0 0.0",
)
transform_map.add_argument(
    "--rot", required=False, default=0.0, type=float, help="rotation in degree"
)
transform_map.add_argument(
    "--axr",
    required=False,
    default=[1, 0, 0],
    nargs="+",
    type=int,
    help="rotation axis",
)
transform_map.add_argument(
    "--out", required=False, default="transformed.mrc", help="output map (mrc/map)"
)

conv_mtz2map = subparsers.add_parser(
    "mtz2map", description="Convert MTZ to MRC/MAP.")
conv_mtz2map.add_argument("--mtz", required=True, help="input map (mtz)")
conv_mtz2map.add_argument("--map", required=True, help="input map (mrc/map)")
conv_mtz2map.add_argument("--out", required=True, help="output map (mrc/map)")

resample_d = subparsers.add_parser(
    "resample", description="Resample a map in Fourier space."
)
resample_d.add_argument("--map", required=True, help="input map (mrc/map)")
resample_d.add_argument(
    "--pix", 
    required=True, 
    nargs="+",
    type=float, 
    help="target pixel size (A)"
)
resample_d.add_argument(
    "--dim",
    required=True,
    #default=None,
    nargs="+",
    type=np.int,
    help="target map dimensions. e.g. 100 100 100 ",
)
resample_d.add_argument(
    "--cel",
    required=False,
    default=None,
    nargs="+",
    type=np.float,
    help="target cell. e.g. a b c 90 90 90",
)
resample_d.add_argument(
    "--out", required=False, default="resampled.mrc", help="output map name"
)

resample_m = subparsers.add_parser(
    "resamplemap2map", description="Resample map2 on map1."
)
resample_m.add_argument("--map1", required=True, help="static map (mrc/map)")
resample_m.add_argument("--map2", required=True,
                        help="map to resample (mrc/map)")
resample_m.add_argument(
    "--out", required=False, default="resampled2staticmap.mrc", help="output map name"
)

realspc = subparsers.add_parser("rcc", description="real space correlation")
realspc.add_argument("--h1", required=True, help="input halfmap1 map")
realspc.add_argument("--h2", required=True, help="input halfmap2 map")
realspc.add_argument("--mdl", required=False, help="Input model (cif/pdb)")
realspc.add_argument("--res", required=False,
                     type=float, help="Resolution (A)")
realspc.add_argument("--msk", required=False, help="input mask (mrc/map)")
realspc.add_argument("--nrm", action="store_true",
                     help="if True use normalized maps")
realspc.add_argument(
    "--knl", required=False, type=int, default=5, help="Kernel size (pixels)"
)
# realspc.add_argument(
#    "--lig", action="store_true", help="use if there is ligand, but no description"
# )
realspc.add_argument(
    "--lgf", required=False, default=None, type=str, help="ligand description file"
)

bfromcc = subparsers.add_parser(
    "bfromcc", description="local b from real space correlation")
bfromcc.add_argument("--h1", required=True, help="input halfmap1 map")
bfromcc.add_argument("--h2", required=True, help="input halfmap2 map")
bfromcc.add_argument("--res", required=True,
                     type=float, help="Resolution (A)")
bfromcc.add_argument("--msk", required=False, help="input mask (mrc/map)")
bfromcc.add_argument(
    "--knl", required=False, type=int, default=5, help="Kernel size (pixels)"
)

mmrealspc = subparsers.add_parser("mmcc", description="real space correlation")
mmrealspc.add_argument("--map", required=True, help="input full/deposited map")
mmrealspc.add_argument("--mdl", required=True, help="Input model (cif/pdb)")
mmrealspc.add_argument("--res", required=True,
                       type=float, help="Resolution (A)")
mmrealspc.add_argument(
    "--nrm", action="store_true", help="if use, normalized maps are used"
)
mmrealspc.add_argument(
    "--msk", required=False, default=None, type=str, help="input mask (mrc/map)"
)
mmrealspc.add_argument(
    "--knl", required=False, type=int, default=5, help="Kernel size (pixels)"
)
mmrealspc.add_argument(
    "--lgf", required=False, default=None, type=str, help="ligand description file"
)

fourierspc = subparsers.add_parser(
    "fcc", description="Fourier space correlation")
fourierspc.add_argument("--h1", required=True, help="input halfmap1 map")
fourierspc.add_argument("--h2", required=True, help="input halfmap2 map")
fourierspc.add_argument(
    "--knl", required=False, type=int, default=5, help="Kernel size (pixels)"
)
fourierspc.add_argument(
    "--msk", required=False, default=None, type=str, help="input mask (mrc/map)"
)

mapmodelvalid = subparsers.add_parser(
    "mapmodelvalidate", description="map-model validation using FSC"
)
mapmodelvalid.add_argument("--h1", required=True, help="input halfmap1 map")
mapmodelvalid.add_argument("--h2", required=True, help="input halfmap2 map")
mapmodelvalid.add_argument("--mdf", required=True,
                           help="input full atomic model")
mapmodelvalid.add_argument(
    "--md1", required=False, default=None, type=str, help="input halfmap1 atomic model"
)
mapmodelvalid.add_argument("--msk", required=False,
                           help="input mask (mrc/map)")
mapmodelvalid.add_argument("--res", required=False,
                           type=float, help="Resolution (A)")
mapmodelvalid.add_argument(
    "--bfc",
    required=False,
    default=0.0,
    type=float,
    help="Overall B factor for model. default=0.0 ",
)
# mapmodelvalid.add_argument(
#    "--lig", action="store_true", help="use if there is ligand, but no description"
# )
mapmodelvalid.add_argument(
    "--lgf", required=False, default=None, type=str, help="ligand description file"
)

mmfsc = subparsers.add_parser("mapmodelfsc", description="map-model FSC")
mmfsc.add_argument("--map", required=True, help="input map")
mmfsc.add_argument("--mdl", required=True, help="input atomic model")
mmfsc.add_argument("--msk", required=False, help="input mask (mrc/map)")
mmfsc.add_argument("--res", required=True, type=float, help="Resolution (A)")
mmfsc.add_argument(
    "--bfc",
    required=False,
    default=0.0,
    type=float,
    help="Overall B factor for model. default=0.0 ignored by REFMAC ",
)
# mmfsc.add_argument(
#    "--lig", action="store_true", help="use if there is ligand, but no description"
# )
mmfsc.add_argument(
    "--lgf", required=False, default=None, type=str, help="ligand description file"
)
mmfsc.add_argument(
    "--phaserand", action="store_true", help="use if phase randomized FSC is calculated"
)

mapoverlay = subparsers.add_parser("overlay", description="overlay maps")
mapoverlay.add_argument(
    "--map", required=True, nargs="+", type=str, help="maplist for overlay"
)
mapoverlay.add_argument(
    "--msk",
    required=False,
    default=None,
    nargs="+",
    type=str,
    help="masklist for overlay",
)
mapoverlay.add_argument(
    "--tra",
    required=False,
    default=[0.0, 0.0, 0.0],
    nargs="+",
    type=float,
    help="translation vector. default=[0.0, 0.0, 0.0]",
)
mapoverlay.add_argument(
    "--rot",
    required=False,
    default=0.0,
    type=float,
    help="rotation in deg. default=0.0",
)
mapoverlay.add_argument(
    "--axr",
    required=False,
    default=[1, 0, 0],
    nargs="+",
    type=int,
    help="rotation axis. default=[1,0,0]",
)
mapoverlay.add_argument(
    "--ncy",
    required=False,
    default=5,
    type=int,
    help="number of fitting cycles. default=5",
)
mapoverlay.add_argument(
    "--res",
    required=False,
    default=6,
    type=float,
    help="starting fit resol. (A). default= 6 A",
)
mapoverlay.add_argument(
    "--fitres",
    required=False,
    default=0.0,
    type=float,
    help="final fit resol. (A). default= 0.0 A",
)
mapoverlay.add_argument(
    "--int",
    required=False,
    default="linear",
    type=str,
    help="interpolation method (linear/cubic). default= linear",
)
mapoverlay.add_argument(
    "--modelres",
    required=False,
    default=5,
    type=float,
    help="model resol. (A). default= 5 A",
)
mapoverlay.add_argument("--hfm", action="store_true",
                        help="if use employ half maps")
mapoverlay.add_argument("--mod", action="store_true",
                        help="if use calls model overlay")
mapoverlay.add_argument("--usecom", action="store_true",
                        help="if used, center-of-mass is used to superimpose maps")


mapaverage = subparsers.add_parser(
    "average", description="weighted average of several maps"
)
mapaverage.add_argument(
    "--map", required=True, nargs="+", type=str, help="maplist to average"
)
mapaverage.add_argument(
    "--msk", required=False, default=None, nargs="+", type=str, help="masklist for maps"
)
mapaverage.add_argument(
    "--tra",
    required=False,
    default=[0.0, 0.0, 0.0],
    nargs="+",
    type=float,
    help="translation vec.",
)
mapaverage.add_argument(
    "--rot", required=False, default=0.0, type=float, help="rotation in deg"
)
mapaverage.add_argument(
    "--axr",
    required=False,
    default=[1, 0, 0],
    nargs="+",
    type=int,
    help="rotation axis",
)
mapaverage.add_argument(
    "--ncy", required=False, default=10, type=int, help="number of fitting cycles"
)
mapaverage.add_argument(
    "--res", required=False, default=6, type=float, help="starting fit resol. (A)"
)
mapaverage.add_argument(
    "--int",
    required=False,
    default="linear",
    type=str,
    help="interpolation method ([linear]/cubic)",
)

diffmap = subparsers.add_parser(
    "diffmap", description="Calculate the difference map"
)
diffmap.add_argument(
    "--map", required=True, nargs="+", type=str, help="maplist to diffmap"
)
diffmap.add_argument(
    "--msk", required=False, default=None, nargs="+", type=str, help="masklist for maps"
)
diffmap.add_argument(
    "--res",
    required=True,
    type=float,
    help="resolution for difference map in Angstroms.",
)
diffmap.add_argument(
    "--mod",
    required=False,
    default="norm",
    type=str,
    help="scaling method. norm (default) - normalized FC, \
        ampli - amplitudes in resolution bins",
)
diffmap.add_argument("--fit", action="store_true",
                     help="if used, maps are superimposed before calculating difference map")
diffmap.add_argument(
    "--ncy", required=False, default=5, type=int, help="number of fitting cycles"
)
diffmap.add_argument(
    "--fitres",
    required=False,
    default=0.0,
    type=float,
    help="final fit resol. (A). default= 0.0 A",
)
diffmap.add_argument("--usecom", action="store_true",
                     help="if used, center-of-mass is used to superimpose maps")

diffmap.add_argument("--usehalfmaps", action="store_true",
                     help="if used, halfmaps are used to calculate difference map")


applymask = subparsers.add_parser(
    "applymask", description="apply mask on the map")
applymask.add_argument("--map", required=True, help="map to be masked")
applymask.add_argument("--msk", required=True, help="mask to be applied")
applymask.add_argument(
    "--out", required=False, default="mapmasked.mrc", help="output map name"
)

scalemap = subparsers.add_parser(
    "scalemap", description="scale onemap to another")
scalemap.add_argument("--m1", required=True, help="input map")
scalemap.add_argument("--m2", required=True, help="map to be scaled")
scalemap.add_argument(
    "--out", required=False, default="scaledmap.mrc", help="output map name"
)

bestmap = subparsers.add_parser("bestmap", description="calculate bestmap")
bestmap.add_argument("--h1", required=True, help="input halfmap1 map")
bestmap.add_argument("--h2", required=True, help="input halfmap2 map")
bestmap.add_argument("--msk", required=False,
                     default=None, help="mask to be applied")
bestmap.add_argument("--B", required=False, type=float, help="relative B")
bestmap.add_argument(
    "--knl", required=False, default=5, type=int, help="kernel radius (pixels)"
)
bestmap.add_argument(
    "--mod",
    required=False,
    default=1,
    type=int,
    help="fsc type (1-resol bins, 2-local)",
)
bestmap.add_argument(
    "--out", required=False, default="bestmap.mrc", help="output map name"
)

predfsc = subparsers.add_parser(
    "predfsc", description="predict FSC based on # particles"
)
predfsc.add_argument("--h1", required=True, help="input halfmap1 map")
predfsc.add_argument("--h2", required=True, help="input halfmap2 map")
predfsc.add_argument("--msk", required=False, help="mask map")
predfsc.add_argument(
    "--npa", required=False, nargs="+", type=float, help="n fold of particles"
)
predfsc.add_argument(
    "--bfc", required=False, nargs="+", type=float, help="list of B factors"
)

refmac = subparsers.add_parser(
    "refmac", description="prepare data for refmac refinement"
)
refmac.add_argument("--h1", required=True, help="input halfmap1 map")
refmac.add_argument("--h2", required=True, help="input halfmap2 map")
refmac.add_argument("--msk", required=False, help="mask map")
refmac.add_argument(
    "--bfc", required=False, nargs="+", type=float, help="b-factor list"
)
refmac.add_argument(
    "--out", required=False, default="output.mtz", help="output mtz file name"
)

occ = subparsers.add_parser(
    "occ", description="overall correlation in real space")
occ.add_argument("--m1", required=True, help="input map1 map")
occ.add_argument("--m2", required=True, help="input map2 map")
occ.add_argument("--msk", required=False, help="mask map")
occ.add_argument(
    "--spc",
    required=False,
    default="real",
    help="space (real/fourier) for CC calculation",
)

mirror = subparsers.add_parser("mirror", description="mirror the map")
mirror.add_argument("--map", required=True, help="input map")

model2map = subparsers.add_parser(
    "model2map", description="calculate model based map")
model2map.add_argument("--mdl", required=True, help="input atomic model")
model2map.add_argument("--res", required=True,
                       type=float, help="Resolution (A)")
model2map.add_argument("--dim", required=True, nargs="+",
                       type=np.int, help="map dim ")
model2map.add_argument(
    "--cel", required=True, nargs="+", type=np.float, help="cell parameters "
)
#model2map.add_argument(
#    "--bfc", required=False, default=0.0, type=float, 
#    help="replace all atomic Bs with this. Only +ve values have an effect"
#)
# model2map.add_argument(
#    "--lig", action="store_true", help="use if there is ligand, but no description"
# )
model2map.add_argument(
    "--lgf", required=False, default=None, type=str, help="ligand description file"
)
model2map.add_argument(
    "--org", required=False, default=None, nargs="+", type=int, help="map origin"
)
model2map.add_argument(
    "--gemmi", action="store_true", 
    help="if used, GEMMI is used instead REFMAC for structure factor calculation"
)

composite = subparsers.add_parser(
    "composite", description="make composite map")
composite.add_argument(
    "--map", required=True, nargs="+", type=str, help="maplist to combine"
)
composite.add_argument(
    "--msk", required=False, default=None, nargs="+", type=str, help="masklist for maps"
)

magref = subparsers.add_parser(
    "magref", description="magnification refinement")
magref.add_argument(
    "--map",
    required=True,
    nargs="+",
    type=str,
    help="maplist to correct for magnification [.mrc/.map]",
)
magref.add_argument("--ref", required=True, type=str,
                    help="reference map [.mrc/.map]")

centerofmass = subparsers.add_parser("com", description="center of mass")
centerofmass.add_argument("--map", required=True,
                          type=str, help="input map (MRC/MAP)")
centerofmass.add_argument(
    "--msk", required=False, default=None, type=str, help="mask to apply on the map"
)

fetchdata = subparsers.add_parser("fetch", description="fetch EMmap and model")
fetchdata.add_argument("--emd", required=True, nargs="+",
                       type=str, help="list of EMD entries. e.g. 3651")
fetchdata.add_argument("--all", action="store_true",
                       help="Use to download all data (mask, map, halfdata, model)")

""" symaxref = subparsers.add_parser(
    "symref", description="refine symmetry axis of a group")
symaxref.add_argument("--map", required=True, nargs="+",
                      type=str, help="list of maps to find symmetry axes")
symaxref.add_argument("--emd", required=False, nargs="+",
                      type=str, help="list of emdbid of maps")
symaxref.add_argument("--res", required=False, nargs="+",
                      type=float, help="list of resolution of maps (A)")
symaxref.add_argument("--mapout", action="store_true",
                      help="Use to output maps") """


pointg = subparsers.add_parser(
    "pointgroup", description="detect point group from the map")
pointg.add_argument("--map", required=True, nargs="+",
                    type=str, help="list of maps to find point groups")
pointg.add_argument("--res", required=True, nargs="+",
                    type=float, help="list of resolution of maps (A)")
pointg.add_argument("--emd", required=False, nargs="+",
                    type=str, help="list of emdbid of maps")
pointg.add_argument("--peak_cutoff", required=False, default=0.8,
                    type=float, help="cutoff for Proshade peak height. default= 0.8")
pointg.add_argument("--use_fsc", action="store_true",
                    help="if used, FSC is used in place for proshade peakheight to decide point group")
pointg.add_argument("--fsc_cutoff", required=False, default=0.7,
                    type=float, help="cutoff for Proshade peak height, default= 0.7")
pointg.add_argument("--ang_tol", required=False, default=5.0,
                    type=float, help="angle tolerence between two axes for determining point group. default= 5 deg.")


symmap = subparsers.add_parser(
    "symmetrise", description="symmetrize map using point group symmetry")
symmap.add_argument("--map", required=True, nargs="+",
                    type=str, help="list of maps to find point groups")
symmap.add_argument("--res", required=True, nargs="+",
                    type=float, help="list of resolution of maps (A)")
symmap.add_argument("--emd", required=False, nargs="+",
                    type=str, help="list of emdbid of maps")
symmap.add_argument("--pointgroup", required=False, nargs="+",
                    type=str, help="list of point groups")
symmap.add_argument("--peak_cutoff", required=False, default=0.8,
                    type=float, help="cutoff for Proshade peak height. default= 0.8")
symmap.add_argument("--use_fsc", action="store_true",
                    help="if used, FSC is used in place for proshade peakheight to decide point group")
symmap.add_argument("--fsc_cutoff", required=False, default=0.7,
                    type=float, help="cutoff for Proshade peak height, default= 0.7")
symmap.add_argument("--ang_tol", required=False, default=5.0,
                    type=float, help="angle tolerence between two axes for determining point group. default= 5 deg.")

rebox = subparsers.add_parser(
    "rebox", description="rebox map and model using a mask")
rebox.add_argument("--map", required=True, nargs="+",
                    type=str, help="list of map names for reboxing")
rebox.add_argument("--msk", required=True, nargs="+",
                    type=str, help="list of mask names for reboxing")
rebox.add_argument("--mdl", required=False, nargs="+",
                    type=str, default=None, help="list of model names (pdb/cif) for reboxing")