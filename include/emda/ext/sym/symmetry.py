"""Map point group detection and symmetry averaging

Author: Rangana Warshamanage
Vesion1: 28.03.2021

This module deals with detecting point group from map and symmetry averaging.

TODO Symmetry averaging in C, D and O groups seem to work fine, but that in
T and I groups may sometime not be optimal as the refined axes can be bit off.

"""

import sys
import numpy as np
import emda.emda_methods as em
from emda.ext.mapfit.utils import double_the_axes
from emda.core.restools import get_resolution_array
from numpy.fft import fftn, ifftn, fftshift, ifftshift
from emda.ext.sym.refine_symaxis import get_pg, refine_pg_generator_axes
from emda.ext.sym.symmetrize_map import symmetrize_map_known_pg, get_matrices, apply_op, rebox_map


class Symmetry:
    """ this class deals with all aspects of map symmetry """
    def __init__(self, imap=None, resol=None, fobj=None):
        self.imap = imap
        self.resol = resol
        self.fobj = fobj
        self.peak_cutoff = 0.8
        self.fsc_cutoff = 0.7
        self.ang_tol = 5.0
        self.fitres = 5.0
        self.fitfsc = 0.7
        self.use_peakheight = True
        self.refine_axes = False
        self.use_fsc = False
        self.emdid = None
        self.point_gp = None
        self.proshade_point_gp = None
        self.generator_orders = None
        self.generator_axes = None
        self.refined_orders = None
        self.refined_axes = None
        self.symavgmap = None
        self.unitcell = None
        self.maporigin = None

    def detect_point_gp(self):
        results = get_pg(
            imap=self.imap,
            resol=self.resol,
            use_peakheight=self.use_peakheight,
            peak_cutoff=self.peak_cutoff,
            use_fsc=self.use_fsc,
            fsc_cutoff=self.fsc_cutoff,
            ang_tol=self.ang_tol,
            fobj=self.fobj,
            )
        if len(results) > 0:
            self.proshade_point_gp = results[0]
            self.point_gp = results[1]
            self.generator_orders = results[2]
            self.generator_axes = results[3]

    def refine_gp_generators(self):
        self.refine_axes = True
        if self.fobj is None:
            self.fobj = open("EMDA_symref.txt", "w")
        results = refine_pg_generator_axes(
            imap=self.imap,
            axlist=self.generator_axes,
            odrlist=self.generator_orders,
            fobj=self.fobj,
            fitres=self.fitres,
            fitfsc=self.fitfsc,
            emdid=self.emdid)
        self.refined_orders = results[0]
        self.refined_axes = results[1]

    def symmetry_average(self):
        if self.point_gp is not None:
            if self.point_gp == 'C1':
                SystemExit("Point group is C1. Can not symmetrize.")
            else:
                results = symmetrize_map_known_pg(
                    imap=self.imap,
                    pg=self.point_gp,
                    outmapname=self.emdid)
        elif self.point_gp is None:
            # detect the point group
            self.detect_point_gp()
            # refine generators
            self.refine_gp_generators()
            # symmetrize map
            if self.point_gp == 'T' or self.point_gp == 'I':
                results = tetrahedral(
                    imap=self.imap,
                    foldlst=self.refined_orders,
                    fnlaxlst=self.refined_axes,
                    outmapname=self.emdid,
                    )
            elif self.point_gp == 'C1':
                SystemExit("Point group is C1. Can not symmetrize.")
            else:
                results = symmetrize_map_known_pg(
                    imap=self.imap,
                    pg=self.point_gp,
                    outmapname=self.emdid)
        self.symavgmap = results[0]
        self.unitcell = results[1]
        self.maporigin = results[2]


def get_pointgroup(maplist, reslist, use_peakheight=True, peak_cutoff=0.8,
                   use_fsc=False, fsc_cutoff=0.7, ang_tol=5.0, emdlist=None,
                   fobj=None):
    if fobj is None:
        fobj = open("EMDA_pointgroup_detection.txt", "w")
    pointgroup_list = []
    proshade_pointgroup_list = []
    for i, imap in enumerate(maplist):
        if i > 0:
            fobj.write("\n")
        obj = Symmetry()
        obj.imap = imap
        obj.fobj = fobj
        if emdlist is not None:
            obj.emdid = emdlist[i]
            fobj.write("EMD-"+str(emdlist[i]))
        obj.resol = reslist[i]
        obj.use_peakheight = use_peakheight
        obj.peak_cutoff = peak_cutoff
        obj.use_fsc = use_fsc
        obj.fsc_cutoff = fsc_cutoff
        obj.ang_tol = ang_tol
        obj.detect_point_gp()
        pointgroup_list.append(obj.point_gp)
        proshade_pointgroup_list.append(obj.proshade_point_gp)
    return pointgroup_list, proshade_pointgroup_list


def tetrahedral(imap, foldlst, fnlaxlst, outmapname=None):
    if outmapname is None:
        outmapname = "emda_sym_averaged_map.mrc"
    ops = get_matrices(foldlst, fnlaxlst)
    uc, arr, orig = em.get_data(imap)
    arr2 = double_the_axes(arr)
    f1 = fftshift(fftn(fftshift(arr2)))
    nbin, _, bin_idx = get_resolution_array(uc, f1)
    frs_sum = f1
    print("Symmetrising map...")
    for op in ops[1:]:
        frs = apply_op(f1, op, bin_idx, nbin)
        frs_sum += frs
    avg_f = frs_sum / len(ops)
    avgmap = ifftshift(np.real(ifftn(ifftshift(avg_f))))
    avgmap = rebox_map(avgmap)
    em.write_mrc(avgmap, outmapname, uc, orig)
    return [avgmap, uc, orig]


def symmetrise_map(maplist, reslist, use_peakheight=True, peak_cutoff=0.8,
                   use_fsc=False, fsc_cutoff=0.7, ang_tol=5.0, pglist=None,
                   emdlist=None, fobj=None):
    if fobj is None:
        fobj = open("EMDA_symmetrisemap.txt", "w")                   
    symavgmaplist = []
    if pglist is not None:
        # do the symmetrisation without pg detection and refinement
        assert len(maplist) == len(reslist) == len(pglist)
        if emdlist is not None:
            assert len(maplist) == len(emdlist)
        for i, imap in enumerate(maplist):
            obj = Symmetry(imap=imap, resol=reslist[i])
            obj.point_gp = pglist[i]
            obj.fobj = fobj
            if emdlist is not None:
                obj.emdid = emdlist[i]
            obj.symmetry_average()
            symavgmaplist.append(obj.symavgmap)
    elif pglist is None:
        # do the symmetrisation after pg detection and refinement
        assert len(maplist) == len(reslist)
        if emdlist is not None:
            assert len(maplist) == len(emdlist)
        for i, imap in enumerate(maplist):
            obj = Symmetry()
            obj.imap = imap
            obj.fobj = fobj
            if emdlist is not None:
                obj.emdid = emdlist[i]
            obj.resol = reslist[i]
            obj.use_peakheight = use_peakheight
            obj.peak_cutoff = peak_cutoff
            obj.use_fsc = use_fsc
            obj.fsc_cutoff = fsc_cutoff
            obj.ang_tol = ang_tol
            obj.symmetry_average()
            symavgmaplist.append(obj.symavgmap)
    return symavgmaplist


if __name__=="__main__":
    map1, resol, pg = sys.argv[1:]

    maplist = []
    maplist.append(map1)
    reslist = []
    reslist.append(float(resol))
    pglist = []
    pglist.append(pg)
    print(maplist, reslist, pglist)
    get_pointgroup(maplist, reslist, use_fsc=True)
    #symmetrise_map(maplist, reslist, pglist=pglist)