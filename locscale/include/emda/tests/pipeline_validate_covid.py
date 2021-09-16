""" Validate COVID structures

This runs several emda tests to help validate
map and map-model fit.
"""

import os, gzip, shutil
import urllib.request
import csv
from contextlib import closing
import datetime
import emda.emda_methods as em
import emda.core.plotter as plotter
import xml.etree.ElementTree as ET


def read_xml(emdid):
    name_list = []
    # read file.xml
    localdb = "/teraraid3/pemsley/emdb/structures/"
    localmmcif = "/net/dstore1/teraraid3/ranganaw/mmCIF/"
    xmlfile1 = "EMD-%s/header/emd-%s.xml" % (emdid, emdid)
    xmlfile2 = "EMD-%s/header/emd-%s-v30.xml" % (emdid, emdid)
    tree1 = ET.parse(localdb + xmlfile1)
    root1 = tree1.getroot()
    tree2 = ET.parse(localdb + xmlfile2)
    root2 = tree2.getroot()

    path = "/net/dstore1/teraraid3/ranganaw/COVID/EMD-%s/" % (emdid)
    try:
        os.mkdir(path)
    except OSError:
        if os.path.isdir(path):
            print("Directory exists!")
        else:
            print("Creation of the directory %s failed" % path)
    else:
        print("Successfully created the directory %s " % path)

    model = None
    maskid = None
    half1id = None
    half2id = None
    claimed_resol = None
    mapname = None
    maskname = None
    half1name = None
    half2name = None
    modelname = None
    for deposition in root1.findall("deposition"):
        for fittedPDBEntryIdList in deposition.findall("fittedPDBEntryIdList"):
            for fittedPDBEntryId in fittedPDBEntryIdList.findall("fittedPDBEntryId"):
                model = fittedPDBEntryId.text
    for map in root1.findall("map"):
        for file in map.findall("file"):
            mapid = file.text
    for processing in root1.findall("processing"):
        for reconstruction in processing.findall("reconstruction"):
            for resolutionByAuthor in reconstruction.findall("resolutionByAuthor"):
                claimed_resol = resolutionByAuthor.text
    for interpretation in root2.findall("interpretation"):
        for segmentation_list in interpretation.findall("segmentation_list"):
            for segmentation in segmentation_list.findall("segmentation"):
                for file in segmentation.findall("file"):
                    maskid = file.text
        for half_map_list in interpretation.findall("half_map_list"):
            for i, half_map in enumerate(half_map_list.findall("half_map")):
                if i == 0:
                    for file in half_map.findall("file"):
                        half1id = file.text
                if i == 1:
                    for file in half_map.findall("file"):
                        half2id = file.text

    print("claimed resol ", claimed_resol)
    print("Mask Id: ", maskid)
    print("Half1id: ", half1id)
    print("Half2id: ", half2id)
    if model is not None:
        """ try:
            cifftplink = "ftp://ftp.ebi.ac.uk/pub/databases/pdb/data/structures/all/mmCIF/{}.cif.gz".format(
                model
            )
            outcif = "%s.cif" % (model)
        except Exception as e:
            print(e)
        try:
            with urllib.request.urlopen(cifftplink) as response:
                with gzip.GzipFile(fileobj=response) as uncompressed:
                    file_content = uncompressed.read()
            with open(path + outcif, "wb") as f:
                f.write(file_content)
        except Exception as e:
            print(e) """
        incif = "%s.cif.gz" % (model)
        outcif = "%s.cif" % (model)
        with gzip.open(localmmcif + incif, "rb") as fmodel:
            file_content = fmodel.read()
        with open(path + outcif, "wb") as f:
            f.write(file_content)        

    readname_list = []
    writename_list = []
    readname_list.append("EMD-%s/map/%s" % (emdid, mapid))
    writename_list.append("emd_%s.map" % (emdid))
    mapname = path + "emd_%s.map" % (emdid)
    if maskid is not None:
        readname_list.append("EMD-%s/masks/%s" % (emdid, maskid))
        writename_list.append("emd_%s_mask.map" % (emdid))
        maskname = path + "emd_%s_mask.map" % (emdid)
    if half1id is not None:
        readname_list.append("EMD-%s/other/%s" % (emdid, half1id))
        writename_list.append("emd_%s_half1.map" % (emdid))
        half1name = path + "emd_%s_half1.map" % (emdid)
        readname_list.append("EMD-%s/other/%s" % (emdid, half2id))
        writename_list.append("emd_%s_half2.map" % (emdid))
        half2name = path + "emd_%s_half2.map" % (emdid)

    for readname, writename in zip(readname_list, writename_list):
        name_list.append(path + writename)
        if readname.endswith((".map")):
            shutil.copy2(localdb + readname, path + writename)
        elif readname.endswith((".gz")):
            with gzip.open(localdb + readname, "rb") as fmap:
                file_content = fmap.read()
            with open(path + writename, "wb") as f:
                f.write(file_content)
    if model is not None:
        name_list.append(path + outcif)
        modelname = path + outcif
    name_list = [mapname, maskname, half1name, half2name, modelname]
    return path, name_list, claimed_resol


def tests(maskname):
    os.chdir(path)
    f = open("EMDA.txt", "w")
    f.write("EMDA session recorded at %s.\n\n" % (datetime.datetime.now()))

    # generate mask if none
    if maskname is None:
        uc, arr, orig = em.get_data(mapname)
        _ = em.mask_from_map(uc=uc, arr=arr, kern=7, itr=2)
        maskname = "mapmask.mrc"

    # generate modelmap if modelname exists
    if modelname is not None:
        uc, arr, orig = em.get_data(mapname)
        modelmap = em.model2map(
            modelxyz=modelname, dim=arr.shape, resol=modelresol, cell=uc, maporigin=orig
        )
        em.write_mrc(modelmap, "modelmap.mrc", unit_cell=uc, map_origin=orig)

    # T1. calculate power spectrum
    if mapname is not None:
        f.write("**** Map power ****\n")
        res_arr, power_spectrum = em.get_map_power(mapname)
        plotter.plot_nlines_log(
            res_arr,
            [power_spectrum],
            curve_label=["Power"],
            mapname="map_power.eps",
            plot_title="Rotationally averaged power spectrum",
        )
        f.write("\n")

    # T2. calculate half map FSC
    if half1name is not None:
        f.write("**** Halfmap FSC ****\n")
        res_arr, fsc_list = em.halfmap_fsc(
            half1name=half1name, half2name=half2name, maskname=maskname
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
        f.write("\n")

    # T3. map-model FSC
    if modelname is not None:
        f.write("**** Map-model FSC ****\n")
        _, _ = em.mapmodel_fsc(
            map1=mapname,
            # model=modelname,
            model="modelmap.mrc",
            mask=maskname,
            modelresol=modelresol,
            fobj=f,
            lig=True,
        )
        f.write("\n")

    # T4. overall cc
    if half1name is not None:
        f.write("**** Overall Correlation ****\n")
        occ, hocc = em.overall_cc(half1name, half2name, space="real", maskname=maskname)
        print("map occ: ", occ, hocc)
        f.write("map occ: " + str(occ))
        f.write("\n")
    if modelname is not None:
        occ, hocc = em.overall_cc(
            map1name=mapname, map2name="modelmap.mrc", space="real", maskname=maskname
        )
        print("map-model occ: ", occ, hocc)
        f.write("map-model occ: " + str(occ))
        f.write("\n")

    # T5. rcc
    if half1name is not None:
        if modelname is not None:
            f.write("**** Local Correlation ****\n")
            em.realsp_correlation(
                half1name,
                half2name,
                kernel_size=kernel_size,
                norm=True,
                # model=modelname,
                model="modelmap.mrc",
                model_resol=modelresol,
                mask_map=maskname,
            )

    # T6. mmcc
    if half1name is None:
        if modelname is not None:
            f.write("**** Map-model Local Correlation ****\n")
            em.realsp_correlation_mapmodel(
                fullmap=mapname,
                # model=modelname,
                model="modelmap.mrc",
                resol=modelresol,
                mask_map=maskname,
                norm=True,
            )

    # T7. bestmap 1dfsc
    if half1name is not None:
        f.write("**** EMDA bestmap ****\n")
        em.bestmap(hf1name=half1name, hf2name=half2name, outfile="bestmap1dfsc.mrc")


if __name__ == "__main__":

    emdbidList = [
        "22251",
        "22659",
        "22532",
        "11639",
        "22083",
        "22292",
        "11719",
        "11526",
        "11728",
        "11617",
        "11328",
        "11329",
        "22293",
        "22821",
        "11144",
        "22256",
        "30276",
        "22221",
        "22078",
        "22156",
        "11068",
        "11119",
        "22729",
        "22127",
        "22161",
        "30374",
        "22491",
        "11497",
        "22352",
        "22301",
        "11174",
        "11205",
        "22162",
        "22158",
        "22253",
        "22254",
        "22255",
        "22256",
        "22660",
        "22668",
        "22533",
        "22534",
        "22535",
        "22574",
        "11728",
        "11526",
        "11330",
        "11331",
        "11332",
        "11333",
        "11334",
        "22822",
        "22823",
        "22824",
        "22825",
        "22826",
        "22831",
        "22832",
        "22833",
        "22834",
        "22835",
        "22836",
        "22837",
        "22838",
        "11145",
        "11146",
        "22251",
        "22253",
        "22254",
        "22255",
        "30277",
        "22222",
        "10863",
        "22730",
        "22731",
        "22732",
        "22733",
        "22734",
        "22735",
        "22736",
        "22737",
        "22124",
        "22125",
        "22126",
        "22128",
        "22492",
        "22494",
        "22497",
        "22506",
        "22507",
        "22508",
        "22512",
        "22516",
        "22517",
        "11493",
        "11494",
        "11495",
        "11496",
        "11498",
        "22353",
        "22354",
        "22355",
        "22356",
        "11173",
        "11184",
        "11203",
        "11204",
        "11206",
        "11207",
        "22159",
        "22275",
    ]
    
    # emdbidList = ["30178"]

    current_path = os.getcwd()

    for emdbid in emdbidList:
        path, name_list, claimed_res = read_xml(emdbid)
        mapname = name_list[0]
        maskname = name_list[1]
        half1name = name_list[2]
        half2name = name_list[3]
        modelname = name_list[4]
        modelresol = float(claimed_res)
        kernel_size = 5
        tests(maskname)
        os.chdir(current_path)
