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
import emda.config
from emda.core import iotools, maptools, restools, plotter, fsc, quaternions

cmdl_parser = argparse.ArgumentParser(
    prog="emda",
    usage="%(prog)s [command] [arguments]",
    formatter_class=argparse.RawDescriptionHelpFormatter,
)

cmdl_parser.add_argument(
    "--version", action="version", version="%(prog)s-" + emda.config.__version__
)

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
######## CCMASK ########
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
######## MAPMASK ########
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
model_mask.add_argument("--binarymask", action="store_true",
                        help="use this to output binary mask")

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
    required=False,
    default=None,
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
    "--modellist",
    required=False,
    default=None,
    nargs="+",
    type=str,
    help="include the list of models on which the \
        transformation (found by map overlay) be applied",
)
mapoverlay.add_argument(
    "--tra",
    required=False,
    default=[0.0, 0.0, 0.0],
    nargs="+",
    type=float,
    help="list of translation vectors. default=[0.0, 0.0, 0.0]",
)
mapoverlay.add_argument(
    "--rot",
    required=False,
    default=[0.0],
    nargs="+",
    type=float,
    help="list of rotations in deg. default=0.0",
)
mapoverlay.add_argument(
    "--axr",
    required=False,
    default=[1, 0, 0],
    nargs="+",
    type=int,
    help="list of rotation axes. default=[1,0,0]",
)
mapoverlay.add_argument(
    "--ncy",
    required=False,
    default=100,
    type=int,
    help="number of fitting cycles. default=5",
)
mapoverlay.add_argument(
    "--fitres",
    required=False,
    default=0.0,
    type=float,
    help="last resolution to use for fitting. default= 0.0 A",
)
mapoverlay.add_argument(
    "--modelres",
    required=False,
    default=None,
    #nargs="+",
    type=float,
    help="resolution for model based map calculation.",
)
mapoverlay.add_argument("--nocom", action="store_true",
                        help="if used, center-of-mass is not used during map overlay")


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
diffmap.add_argument("--nocom", action="store_true",
                     help="if used, center-of-mass is not used during overlay")

diffmap.add_argument("--usehalfmaps", action="store_true",
                     help="if used, halfmaps are used to calculate difference map")
diffmap.add_argument(
    "--tra",
    required=False,
    default=[0.0, 0.0, 0.0],
    nargs="+",
    type=float,
    help="list of translation vectors. default=[0.0, 0.0, 0.0]",
)
diffmap.add_argument(
    "--rot",
    required=False,
    default=[0.0],
    nargs="+",
    type=float,
    help="list of rotations in deg. default=0.0",
)
diffmap.add_argument(
    "--axr",
    required=False,
    default=[1, 0, 0],
    nargs="+",
    type=int,
    help="list of rotation axes. default=[1,0,0]",
)


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
model2map.add_argument(
    "--shift_to_boxcenter", action="store_true", 
    help="if used, calculated map is placed at the boxcenter"
)

composite = subparsers.add_parser(
    "composite", description="make composite map")
composite.add_argument(
    "--map", required=True, nargs="+", type=str, help="maplist to combine"
)
composite.add_argument(
    "--msk", required=False, default=None, nargs="+", type=str, help="masklist for maps"
)

magref = subparsers.add_parser("magref", description="magnification refinement")
magref.add_argument(
                    "--map",
                    required=True,
                    nargs="+",
                    type=str,
                    help="maplist to correct for magnification [.mrc/.map]")
magref.add_argument("--ref", 
                    required=True, 
                    type=str,
                    help="reference map [.mrc/.map]")
magref.add_argument("--res", 
                    required=False, 
                    type=float,
                    default=4.0,
                    help="data resolution for magnification refinement.\
                        This limit will be imposed on all maps.")
magref.add_argument("--msk",
                    required=False,
                    default=None,
                    nargs="+",
                    type=str,
                    help="list of masks to apply on maps[.mrc/.map]")

centerofmass = subparsers.add_parser("com", description="center of mass")
centerofmass.add_argument("--map", 
                          required=True,
                          type=str, help="input map (MRC/MAP)")
centerofmass.add_argument(
                          "--msk", 
                          required=False, 
                          default=None, 
                          type=str, 
                          help="mask to apply on the map")

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


def apply_mask(args):
    from emda.emda_methods import applymask

    applymask(args.map, args.msk, args.out)


def map_info(args):
    from emda.emda_methods import read_map

    uc, arr, origin = read_map(args.map)
    print("Unit cell: ", uc)
    print("Sampling: ", arr.shape)
    print("Pixel size: ", round(uc[0] / arr.shape[0], 3))
    print("Origin: ", origin)


def anymap_fsc(args, fobj):
    from emda.emda_methods import twomap_fsc

    res_arr, bin_fsc = twomap_fsc(args.map1, args.map2, fobj=fobj)
    plotter.plot_nlines(
        res_arr, [bin_fsc], "twomap_fsc.eps", curve_label=["twomap_fsc"]
    )


def halfmap_fsc(args):
    from emda.emda_methods import halfmap_fsc, halfmap_fsc_ph

    if args.phaserand:
        res_arr, fsc_list = halfmap_fsc_ph(
            half1name=args.h1, half2name=args.h2, filename=args.out, maskname=args.msk)
    else:
        res_arr, fsc_list = halfmap_fsc(
            half1name=args.h1, half2name=args.h2, filename=args.out, maskname=args.msk
        )
        if len(fsc_list) == 2:
            plotter.plot_nlines(
                res_arr,
                fsc_list,
                "halfmap_fsc.eps",
                curve_label=["unmask-FSC", "masked-FSC"],
                plot_title="Halfmap FSC",
            )
        elif len(fsc_list) == 1:
            plotter.plot_nlines(
                res_arr,
                fsc_list,
                "halfmap_fsc.eps",
                curve_label=["unmask-FSC"],
                plot_title="Halfmap FSC",
            )


def singlemap_fsc(args):
    from emda.emda_methods import singlemap_fsc as sfsc

    res_arr, bin_fsc, _ = sfsc(map1name=args.h1, knl=args.knl)
    plotter.plot_nlines(res_arr, [bin_fsc],
                        "map_fsc.eps", curve_label=["map_fsc"])


def cc_mask(args):
    from emda import emda_methods as em

    maskname = "halfmap_mask.mrc"
    uc, arr1, origin = em.read_map(args.h1)
    uc, arr2, origin = em.read_map(args.h2)
    ccmask = em.mask_from_halfmaps(
        uc=uc,
        half1=arr1,
        half2=arr2,
        radius=args.knl,
        iter=args.itr,
        dthresh=args.dthreshold,
    )
    em.write_mrc(ccmask, maskname, uc, origin)


def lowpass_map(args):
    from emda import emda_methods as em

    uc, map1, orig = em.read_map(args.map)
    _, map_lwp = em.lowpass_map(
        uc=uc, arr1=map1, resol=args.res, filter=args.fil)
    if args.fil == "butterworth":
        outname = "{0}_{1}.{2}".format("lowpass_bw", str(args.res), "mrc")
    else:
        outname = "{0}_{1}.{2}".format("lowpass", str(args.res), "mrc")
    em.write_mrc(map_lwp, outname, uc, orig)


def power_map(args):
    from emda.emda_methods import get_map_power

    res_arr, power_spectrum = get_map_power(args.map)
    plotter.plot_nlines_log(
        res_arr,
        [power_spectrum],
        curve_label=["Power"],
        mapname="map_power.eps",
        plot_title="Rotationally averaged power spectrum",
    )


def mapresol(args):
    from emda.emda_methods import estimate_map_resol

    resol = estimate_map_resol(args.h1, args.h2)
    print("Map resolution (A):", resol)


def map2mtz(args):
    from emda import emda_methods

    if args.out.endswith((".mtz")):
        outfile = args.out
    else:
        outfile = args.out + ".mtz"
    emda_methods.map2mtz(args.map, outfile, resol=args.res)


def mtz2map(args):
    from emda import emda_methods

    uc, ar, origin = iotools.read_map(args.map)
    dat = emda_methods.mtz2map(args.mtz, ar.shape)
    if args.out.endswith((".mrc")):
        outfile = args.out
    else:
        outfile = args.out + ".mrc"
    iotools.write_mrc(dat, outfile, uc, origin)


#def resample_data(args):
#    import numpy as np
#    from emda.emda_methods import read_map, resample_data, write_mrc
#    import emda.ext.mapfit.utils as utils
#
#    ucx, arr, org = read_map(args.map)
#    uc = np.zeros(3, 'float')
#    uc[0] = ucx[2]
#    uc[1] = ucx[1]
#    uc[2] = ucx[0]
#    arr = utils.set_dim_even(arr)
#    arr = np.transpose(arr)
#    curnt_pix = float(round(uc[0] / arr.shape[0], 5))
#    if args.cel:
#        target_uc = args.cel
#    if args.pix is None:
#        targt_pix = curnt_pix
#    else:
#        targt_pix = float(round(args.pix, 5))
#    # print('pixel size [current, target]: ', curnt_pix,targt_pix)
#    if args.dim is None:
#        if args.cel:
#            dim1 = int(round(target_uc[0] / targt_pix))
#            dim2 = int(round(target_uc[1] / targt_pix))
#            dim3 = int(round(target_uc[2] / targt_pix))
#            new_arr = resample_data(
#                curnt_pix=curnt_pix,
#                targt_pix=targt_pix,
#                targt_dim=[dim1, dim2, dim3],
#                arr=arr,
#            )
#        else:
#            dim1 = int(round(uc[0] / targt_pix))
#            dim2 = int(round(uc[1] / targt_pix))
#            dim3 = int(round(uc[2] / targt_pix))
#            new_arr = resample_data(
#                curnt_pix=curnt_pix,
#                targt_pix=targt_pix,
#                targt_dim=[dim1, dim2, dim3],
#                arr=arr,
#            )
#            new_arr = np.transpose(new_arr)
#            target_uc = round(targt_pix, 3) * np.asarray([dim1, dim2, dim3], dtype="int")
#            #target_uc = uc
#    if args.dim is not None:
#        if args.cel:
#            if abs(targt_pix - round(target_uc[0] / args.dim[0], 3)) < 10e-3:
#                new_arr = resample_data(
#                    curnt_pix=curnt_pix,
#                    targt_pix=targt_pix,
#                    targt_dim=args.dim,
#                    arr=arr,
#                )
#            else:
#                print(
#                    "target pixel size does not match \
#                    with given cell and dims."
#                )
#                exit()
#        else:
#            target_uc = round(targt_pix, 3) * np.asarray(args.dim, dtype="int")
#            print("New cell: ", target_uc)
#            new_arr = resample_data(
#                curnt_pix=curnt_pix, targt_pix=targt_pix, targt_dim=args.dim, arr=arr
#            )
#    write_mrc(new_arr, args.out, target_uc, org)


def resample_data(args):
    import numpy as np
    from emda.emda_methods import read_map, resample_data, write_mrc
    import emda.ext.mapfit.utils as utils

    uc, arr, org = read_map(args.map)
    #arr = utils.set_dim_even(arr)
    cpix1 = float(round(uc[0] / arr.shape[0], 5))
    cpix2 = float(round(uc[1] / arr.shape[1], 5))
    cpix3 = float(round(uc[2] / arr.shape[2], 5))
    curnt_pix = [cpix1, cpix2, cpix3]
    if args.pix is None:
        targt_pix = curnt_pix
    else:
        targt_pix = args.pix
    new_arr = resample_data(
        curnt_pix=curnt_pix,
        targt_pix=targt_pix,
        targt_dim=args.dim,
        arr=arr,
    )
    target_uc = np.array(targt_pix, 'float') * np.asarray(new_arr.shape, dtype="int")
    write_mrc(new_arr, args.out, target_uc, org)


def resample2maps(args):
    # Resampling map2 on map1
    from emda.emda_methods import read_map, resample_data, write_mrc
    import emda.ext.mapfit.utils as utils

    uc1, arr1, org1 = read_map(args.map1)
    uc2, arr2, org2 = read_map(args.map2)
    arr1 = utils.set_dim_even(arr1)
    arr2 = utils.set_dim_even(arr2)
    curnt_pix, targt_pix = [], []
    for i in range(3):
        curnt_pix.append(float(round(uc2[0] / arr2.shape[0], 5))) 
        targt_pix.append(float(round(uc1[0] / arr1.shape[0], 5))) 
    new_arr = resample_data(
        curnt_pix=curnt_pix, targt_pix=targt_pix, targt_dim=arr1.shape, arr=arr2
    )
    write_mrc(new_arr, args.out, uc1, org1)


def realsp_corr(args):
    from emda.emda_methods import realsp_correlation

    realsp_correlation(
        half1map=args.h1,
        half2map=args.h2,
        kernel_size=args.knl,
        norm=args.nrm,
        model=args.mdl,
        model_resol=args.res,
        mask_map=args.msk,
        lgf=args.lgf,
    )


def b_from_cc(args):
    from emda.emda_methods import b_from_correlation

    b_from_correlation(
        half1map=args.h1,
        half2map=args.h2,
        kernel_size=args.knl,
        resol=args.res,
        mask_map=args.msk,
    )


def mmrealsp_corr(args):
    from emda.emda_methods import realsp_correlation_mapmodel

    realsp_correlation_mapmodel(
        fullmap=args.map,
        kernel_size=args.knl,
        model=args.mdl,
        resol=args.res,
        mask_map=args.msk,
        norm=args.nrm,
        lgf=args.lgf,
    )


def fouriersp_corr(args):
    from emda.emda_methods import fouriersp_correlation

    # half1_map, half2_map, kernel_size=5, mask=None
    fouriersp_correlation(
        half1_map=args.h1, half2_map=args.h2, kernel_size=args.knl, mask=args.msk
    )


def validate_mapmodel(args):
    from emda.emda_methods import map_model_validate

    _ = map_model_validate(
        half1map=args.h1,
        half2map=args.h2,
        modelfpdb=args.mdf,
        model1pdb=args.md1,
        mask=args.msk,
        modelresol=args.res,
        bfac=args.bfc,
        # lig=args.lig,
        lgf=args.lgf,
    )


def mapmodel_fsc(args, fobj):
    import emda.emda_methods as em

    _, _ = em.mapmodel_fsc(
        map1=args.map,
        model=args.mdl,
        bfac=args.bfc,
        # lig=args.lig,
        mask=args.msk,
        modelresol=args.res,
        lgf=args.lgf,
        phaserand=args.phaserand,
        fobj=fobj,
    )


def map_overlay(args):
    from emda.emda_methods import overlay_maps

    overlay_maps(
        maplist=args.map,
        masklist=args.msk,
        modellist=args.modellist,
        tra=args.tra,
        rot=args.rot,
        axr=args.axr,
        ncy=args.ncy,
        modelres=args.modelres,
        fitres=args.fitres,
        nocom=args.nocom,
    )


def map_transform(args):
    from emda.emda_methods import map_transform

    map_transform(args.map, args.tra, args.rot, args.axr, args.out)


def map_average(args, fobj):
    from emda.emda_methods import average_maps

    fobj.write("***** Map Average *****\n")
    average_maps(
        maplist=args.map,
        masklist=args.msk,
        tra=args.tra,
        rot=args.rot,
        axr=args.axr,
        ncy=args.ncy,
        res=args.res,
        fobj=fobj,
        interp=args.int,
    )


def apply_bfac(args):
    from emda.emda_methods import apply_bfactor_to_map

    all_maps = apply_bfactor_to_map(args.map, args.bfc, args.out)


def half_to_full(args):
    from emda.emda_methods import half2full

    fullmap = half2full(args.h1, args.h2, args.out)


def diff_map(args):
    from emda.emda_methods import difference_map

    difference_map(maplist=args.map, masklist=args.msk,
                   diffmapres=args.res, mode=args.mod, ncy=args.ncy,
                   fit=args.fit, usehalfmaps=args.usehalfmaps, 
                   usecom=args.nocom, fitres=args.fitres, rot=args.rot,
                   axr=args.axr, trans=args.tra, )


def scale_map(args):
    from emda.emda_methods import scale_map2map

    scaled_map = scale_map2map(args.m1, args.m2, args.out)


def best_map(args):
    from emda.emda_methods import bestmap

    bestmap(
        hf1name=args.h1,
        hf2name=args.h2,
        outfile=args.out,
        mode=args.mod,
        knl=args.knl,
        mask=args.msk,
        B=args.B,
    )


def pred_fsc(args):
    from emda.emda_methods import predict_fsc
    import numpy as np

    # npa = np.asarray(args.npa, dtype='float')
    fsc_lst, res_arr, bin_idx, nbin = predict_fsc(
        hf1name=args.h1,
        hf2name=args.h2,
        nparticles=args.npa,
        bfac=args.bfc,
        mask=args.msk,
    )


def refmac_data(args):
    from emda.emda_methods import prepare_refmac_data

    prepare_refmac_data(
        hf1name=args.h1,
        hf2name=args.h2,
        maskname=args.msk,
        bfac=args.bfc,
        outfile=args.out,
    )


def maptomtzfull(args):
    from emda.emda_methods import read_map, map2mtzfull

    uc, arr1, orig = read_map(args.h1)
    uc, arr2, orig = read_map(args.h2)
    map2mtzfull(uc=uc, arr1=arr1, arr2=arr2, mtzname=args.out)


def overallcc(args):
    from emda.emda_methods import overall_cc

    occ, hocc = overall_cc(
        map1name=args.m1, map2name=args.m2, maskname=args.msk, space=args.spc
    )


def mirrormap(args):
    from emda.emda_methods import mirror_map

    mirror_map(args.map)


def modeltomap(args):
    from emda.emda_methods import model2map, model2map_gm, write_mrc

    if args.gemmi:
        # GEMMI for sfcalc
        modelmap = model2map_gm(
            modelxyz=args.mdl, 
            resol=args.res,
            dim=args.dim, 
            cell=args.cel, 
            maporigin=args.org,
            shift_to_boxcenter=args.shift_to_boxcenter,
            )
        write_mrc(modelmap, "modelmap_gm.mrc", args.cel, args.org)
    else:
        # REFMAC sfcalc
        modelmap = model2map(
            modelxyz=args.mdl,
            dim=args.dim,
            resol=args.res,
            cell=args.cel,
            maporigin=args.org,
            # lig=args.lig,
            ligfile=args.lgf,
            shift_to_boxcenter=args.shift_to_boxcenter,
        )
        write_mrc(modelmap, "modelmap_refmac.mrc", args.cel, args.org)
    


def mask4mmap(args):
    from emda import emda_methods as em

    uc, arr, orig = em.read_map(args.map)
    _ = em.mask_from_map(
        uc=uc,
        arr=arr,
        orig=orig,
        kern=args.knl,
        resol=args.res,
        filter=args.fil,
        prob=args.prb,
        itr=args.itr,
    )


def mask4mmodel(args):
    from emda import emda_methods as em

    _ = em.mask_from_atomic_model(
        mapname=args.map, modelname=args.mdl, atmrad=args.atmrad, 
        binary_mask=args.binarymask)


def composite_map(args):
    from emda import emda_methods as em

    em.compositemap(maps=args.map, masks=args.msk)


def magnification(args):
    from emda import emda_methods as em

    em.mapmagnification(maplist=args.map, 
                        rmap=args.ref, 
                        masklist=args.msk, 
                        resol=args.res)


def center_of_mass(args):
    from emda import emda_methods as em

    uc, arr, orig = em.get_data(args.map)
    if args.msk is not None:
        _, mask, _ = em.get_data(args.msk)
        assert arr.shape == mask.shape
        arr = arr * mask
    print(em.center_of_mass_density(arr))


def fetch_data(args):
    from emda import emda_methods as em
    em.fetch_data(args.emd, args.all)


""" def symaxis_refinement(args):
    from emda import emda_methods as em

    _, _, _, _, _ = em.symaxis_refine(
        maplist=args.map, emdbidlist=args.emd, mapoutvar=args.mapout, reslist=args.res) """


def pointgroup(args, fobj):
    from emda import emda_methods as em

    _, _ = em.get_map_pointgroup(
        maplist=args.map,
        reslist=args.res,
        use_peakheight=True,  # args.use_peakheight,
        peak_cutoff=args.peak_cutoff,
        use_fsc=args.use_fsc,
        fsc_cutoff=args.fsc_cutoff,
        ang_tol=args.ang_tol,
        emdlist=args.emd,
        fobj=fobj,
    )


def symmetrize(args, fobj):
    from emda import emda_methods as em

    _ = em.symmetry_average(
        maplist=args.map,
        reslist=args.res,
        use_peakheight=True,  # args.use_peakheight,
        peak_cutoff=args.peak_cutoff,
        use_fsc=args.use_fsc,
        fsc_cutoff=args.fsc_cutoff,
        ang_tol=args.ang_tol,
        emdlist=args.emd,
        pglist=args.pointgroup,
        fobj=fobj,
    )


def reboxmap(args):
    from emda import emda_methods as em

    em.rebox_mapmodel(maplist=args.map, masklist=args.msk, modellist=args.mdl)


def main(command_line=None):
    f = open("EMDA.txt", "w")
    f.write("EMDA session recorded at %s.\n\n" % (datetime.datetime.now()))
    args = cmdl_parser.parse_args(command_line)
    if args.command == "info":
        map_info(args)
    if args.command == "fsc":
        anymap_fsc(args, f)
        f.close()
    if args.command == "halffsc":
        halfmap_fsc(args)
    if args.command == "ccmask":
        cc_mask(args)
    if args.command == "lowpass":
        lowpass_map(args)
    if args.command == "power":
        power_map(args)
    if args.command == "resol":
        mapresol(args)
    if args.command == "map2mtz":
        map2mtz(args)
    if args.command == "mtz2map":
        mtz2map(args)
    if args.command == "resample":
        resample_data(args)
    if args.command == "rcc":
        realsp_corr(args)
    if args.command == "bfromcc":
        b_from_cc(args)
    if args.command == "mmcc":
        mmrealsp_corr(args)
    if args.command == "fcc":
        fouriersp_corr(args)
    if args.command == "mapmodelvalidate":
        validate_mapmodel(args)
    if args.command == "mapmodelfsc":
        mapmodel_fsc(args, f)
        f.close()
    if args.command == "overlay":
        map_overlay(args)
    if args.command == "average":
        map_average(args, f)
        f.close()
    if args.command == "transform":
        map_transform(args)
    if args.command == "bfac":
        apply_bfac(args)
    if args.command == "singlemapfsc":
        singlemap_fsc(args)
    if args.command == "half2full":
        half_to_full(args)
    if args.command == "diffmap":
        diff_map(args)
    if args.command == "applymask":
        apply_mask(args)
    if args.command == "scalemap":
        scale_map(args)
    if args.command == "bestmap":
        best_map(args)
    if args.command == "predfsc":
        pred_fsc(args)
    if args.command == "refmac":
        refmac_data(args)
    if args.command == "occ":
        overallcc(args)
    if args.command == "mirror":
        mirrormap(args)
    if args.command == "model2map":
        modeltomap(args)
    if args.command == "map2mtzfull":
        maptomtzfull(args)
    if args.command == "mapmask":
        mask4mmap(args)
    if args.command == "modelmask":
        mask4mmodel(args)
    if args.command == "composite":
        composite_map(args)
    if args.command == "resamplemap2map":
        resample2maps(args)
    if args.command == "magref":
        magnification(args)
    if args.command == "com":
        center_of_mass(args)
    if args.command == "fetch":
        fetch_data(args)
    """ if args.command == "symref":
        symaxis_refinement(args) """
    if args.command == "pointgroup":
        pointgroup(args, f)
    if args.command == "symmetrise":
        symmetrize(args, f)
    if args.command == "rebox":
        reboxmap(args)


if __name__ == "__main__":
    main()
