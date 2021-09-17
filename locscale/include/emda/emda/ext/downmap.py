import os, gzip, shutil
import urllib.request
from contextlib import closing
import xml.etree.ElementTree as ET


def fetch_data(emdbid):
    headerxml = (
        "ftp://ftp.ebi.ac.uk/pub/databases/emdb/structures/EMD-%s/header/emd-%s.xml"
        % (emdbid, emdbid)
    )
    header30xml = (
        "ftp://ftp.ebi.ac.uk/pub/databases/emdb/structures/EMD-%s/header/emd%s-v30.xml"
        % (emdbid, emdbid)
    )

    outmap = None; outcif = None; claimed_resol = None

    with closing(urllib.request.urlopen(headerxml)) as r:
        with open("file.xml", "wb") as f:
            shutil.copyfileobj(r, f)

    tree = ET.parse("file.xml")
    root = tree.getroot()

    model = None
    for deposition in root.findall("deposition"):
        for fittedPDBEntryIdList in deposition.findall("fittedPDBEntryIdList"):
            for fittedPDBEntryId in fittedPDBEntryIdList.findall("fittedPDBEntryId"):
                model = fittedPDBEntryId.text

    for processing in root.findall("processing"):
        for reconstruction in processing.findall("reconstruction"):
            for resolutionByAuthor in reconstruction.findall("resolutionByAuthor"):
                claimed_resol = resolutionByAuthor.text

    # create a directory with EMDBID
    #path = "/Users/sandy/MRC/REFMAC/COVID19/EMD-%s/" % (emdbid)
    path = os.getcwd() + "/EMD-%s/" % (emdbid)
    try:
        os.mkdir(path)
    except OSError:
        if os.path.isdir(path):
            print("Directory exists!")
        else:
            print("Creation of the directory %s failed" % path)
    else:
        print("Successfully created the directory %s " % path)

    if model is not None:
        try:
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
            print(e)
    mapftplink = (
        "ftp://ftp.ebi.ac.uk/pub/databases/emdb/structures/EMD-%s/map/emd_%s.map.gz"
        % (emdbid, emdbid)
    )

    outmap = "emd_%s.map" % (emdbid)
    try:
        with urllib.request.urlopen(mapftplink) as response:
            with gzip.GzipFile(fileobj=response) as uncompressed:
                file_content = uncompressed.read()
        with open(path + outmap, "wb") as f:
            f.write(file_content)
    except Exception as e:
        print(e)


def fetch_all_data(emdid):
    name_list = []
    headerxml = (
        "ftp://ftp.ebi.ac.uk/pub/databases/emdb/structures/EMD-%s/header/emd-%s.xml"
        % (emdid, emdid)
    )
    header30xml = (
        "ftp://ftp.ebi.ac.uk/pub/databases/emdb/structures/EMD-%s/header/emd-%s-v30.xml"
        % (emdid, emdid)
    )
    with closing(urllib.request.urlopen(headerxml)) as r:
        with open("file.xml", "wb") as f:
            shutil.copyfileobj(r, f)
    tree1 = ET.parse("file.xml")
    root1 = tree1.getroot()

    with closing(urllib.request.urlopen(header30xml)) as r:
        with open("file30.xml", "wb") as f:
            shutil.copyfileobj(r, f)
    tree2 = ET.parse("file30.xml")
    root2 = tree2.getroot()
    
    path = os.getcwd() + "/EMD-%s/" % (emdid)
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
        try:
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
            print(e)

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

    ftplink = (
        "ftp://ftp.ebi.ac.uk/pub/databases/emdb/structures/"
    )
    for readname, writename in zip(readname_list, writename_list):
        name_list.append(path + writename)
        try:
            fid = ftplink+readname
            if fid.endswith((".map")):
                with closing(urllib.request.urlopen(fid)) as r:
                    with open(path + writename, 'wb') as f:
                        shutil.copyfileobj(r, f)
            elif fid.endswith((".gz")):
                with urllib.request.urlopen(ftplink+readname) as response:
                    with gzip.GzipFile(fileobj=response) as uncompressed:
                        file_content = uncompressed.read()
                with open(path + writename, "wb") as f:
                    f.write(file_content)
        except Exception as e:
            print(e)


def main(emdbidList, alldata=False):
    for emdbid in emdbidList:
        if alldata:
            fetch_all_data(emdbid)
        else:
            fetch_data(emdbid)
        print("Fetched ", emdbid)

""" if __name__ == "__main__":
    emdbidList = ["11310"]#, "11276", "11288", "11321", "11325"]
    for emdbid in emdbidList:
        fetch_data(emdbid)
    print("Fetched data") """