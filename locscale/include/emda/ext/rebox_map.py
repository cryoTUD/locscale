"""Rebox map, mask and model"""

import numpy as np
import emda.emda_methods as em


def get_minimum_bounding_dims(arr1):
    from scipy import ndimage

    # box dim check begin
    indices = np.indices(arr1.shape)
    map_sum = np.sum(arr1) / 1000
    print(map_sum)
    map_sum = 0.
    arr1 = arr1 * (arr1 > map_sum)
    com = ndimage.measurements.center_of_mass(arr1)
    dist2 = (
          (indices[0, :, :] - com[0]) ** 2
        + (indices[1, :, :] - com[1]) ** 2
        + (indices[2, :, :] - com[2]) ** 2
    )
    # np.max(dist2[mask > 0.5])
    min_dim = int(np.sqrt(np.max(dist2)) * 2) + 1
    print(com, min_dim)
    # box dim check end
    return com, min_dim


def get_reboxed_cubic_image(arr, mask):
    #arr = arr * mask
    mask = mask * (mask > 1.e-5)
    i, j, k = np.nonzero(mask)
    z2, y2, x2 = np.max(i), np.max(j), np.max(k)
    z1, y1, x1 = np.min(i), np.min(j), np.min(k)
    dimz = z2 - z1
    dimy = y2 - y1
    dimx = x2 - x1
    dim = np.max([dimz, dimy, dimx])
    if dim % 2 != 0:
        dim += 1
    newarr  = np.zeros((dim+10, dim+10, dim+10), 'float')
    newmask = np.zeros((dim+10, dim+10, dim+10), 'float')
    print((dim - dimz), (dim - dimy), (dim - dimx))
    dz = (dim - dimz) // 2
    dy = (dim - dimy) // 2
    dx = (dim - dimx) // 2
    dz += 5
    dy += 5
    dx += 5
    newarr[dz:dz+dimz, dy:dy+dimy, dx:dx+dimx] = arr[z1:z2, y1:y2, x1:x2]
    newmask[dz:dz+dimz, dy:dy+dimy, dx:dx+dimx] = mask[z1:z2, y1:y2, x1:x2]
    return newarr, newmask


def model_rebox(arr, mask, mmcif_file, uc=None):
    import gemmi
    import os
    from emda.core.iotools import pdb2mmcif

    arr = arr * mask
    mask = mask * (mask > 1.e-5)
    i, j, k = np.nonzero(mask)
    z2, y2, x2 = np.max(i), np.max(j), np.max(k)
    z1, y1, x1 = np.min(i), np.min(j), np.min(k)
    dimz = z2 - z1
    dimy = y2 - y1
    dimx = x2 - x1
    if mmcif_file.endswith(".pdb"):
        pdb2mmcif(mmcif_file)
        mmcif_file = "./out.cif"
    doc = gemmi.cif.read_file(mmcif_file)
    block = doc.sole_block()  # cif file as a single block
    col_x = block.find_values("_atom_site.Cartn_x")
    col_y = block.find_values("_atom_site.Cartn_y")
    col_z = block.find_values("_atom_site.Cartn_z")
    if uc is not None:
        a, b, c = uc[2], uc[1], uc[0]
        alf = bet = gam = 90.
    else:
        a = block.find_value("_cell.length_a")
        b = block.find_value("_cell.length_b")
        c = block.find_value("_cell.length_c")
        alf = block.find_value("_cell.angle_alpha")
        bet = block.find_value("_cell.angle_beta")
        gam = block.find_value("_cell.angle_gamma")
    cell = np.array([a, b, c, alf, bet, gam], dtype="float")
    pixc = float(c) / arr.shape[0]
    pixb = float(b) / arr.shape[1]
    pixa = float(a) / arr.shape[2]
    cart_z = pixc * z1
    cart_y = pixb * y1
    cart_x = pixa * x1
    dim = np.max([dimz, dimy, dimx])
    if dim % 2 != 0:
        dim += 1
    dz = (dim - dimz) // 2
    dy = (dim - dimy) // 2
    dx = (dim - dimx) // 2
    for n, _ in enumerate(col_x):
        col_x[n] = str((float(col_x[n]) - cart_x) + (dx + 5) * pixa)
        col_y[n] = str((float(col_y[n]) - cart_y) + (dy + 5) * pixb)
        col_z[n] = str((float(col_z[n]) - cart_z) + (dz + 5) * pixc)
    doc.write_file("./tmp.cif")
    st = gemmi.read_structure("./tmp.cif")
    ca = (dim+10) * pixa
    cb = (dim+10) * pixb
    cc = (dim+10) * pixc
    st.cell.set(ca, cb, cc, 90., 90., 90.)
    st.make_mmcif_document().write_file("emda_reboxed_model.cif")
    os.remove("./out.cif")
    os.remove("./tmp.cif")


def reboxmap(arr, uc, mask, outname="reboxedmap.mrc"):
    pixs1 = uc[0]/arr.shape[0]
    pixs2 = uc[1]/arr.shape[1]
    pixs3 = uc[2]/arr.shape[2]
    newarr, newmask = get_reboxed_cubic_image(arr, mask)
    dim = newarr.shape[0]
    uc2 = [pixs1*dim, pixs2*dim, pixs3*dim]
    em.write_mrc(newarr, outname, uc2)
    em.write_mrc(newmask, "emda_reboxedmask.mrc", uc2)
    return newarr, uc2


def _rebox_using_maponly(imap):
    uc, arr, orig = em.get_data(imap)
    pixs = uc[0]/arr.shape[0]
    _, arrlp = em.lowpass_map(uc, arr, 15., 'butterworth', order=1)
    _, dim = get_minimum_bounding_dims(arrlp)
    com = em.center_of_mass_density(arr)
    nz, ny, nx = arr.shape
    box_center = (nz//2, ny//2, nx//2)
    shift = np.subtract(com, box_center)
    arr_mvd = em.shift_density(arr, shift)
    dx = (nx - dim)//2
    uc2 = [pixs*dim, pixs*dim, pixs*dim]
    reboxedmap = arr[dx:dx+dim, dx:dx+dim, dx:dx+dim]
    em.write_mrc(reboxedmap, "reboxedmap_xx.mrc", uc2)
    return reboxedmap


def rebox_model(mmcif_file, uc_old, uc_new):
    import gemmi
    import numpy as np
    
    doc = gemmi.cif.read_file(mmcif_file)
    st = gemmi.read_structure(mmcif_file)
    model = st[0]
    com = model.calculate_center_of_mass()
    #print(com)
    block = doc.sole_block()  # cif file as a single block
    col_x = block.find_values("_atom_site.Cartn_x")
    col_y = block.find_values("_atom_site.Cartn_y")
    col_z = block.find_values("_atom_site.Cartn_z")
    a = block.find_value("_cell.length_a")
    b = block.find_value("_cell.length_b")
    c = block.find_value("_cell.length_c")
    alf = block.find_value("_cell.angle_alpha")
    bet = block.find_value("_cell.angle_beta")
    gam = block.find_value("_cell.angle_gamma")
    cell = np.array([a, b, c, alf, bet, gam], dtype="float")
    pixa = uc_new[0] / uc_old[0]
    pixb = uc_new[1] / uc_old[1]
    pixc = uc_new[2] / uc_old[2]
    for n, _ in enumerate(col_x):
        col_x[n] = str(float(col_x[n]) * pixa)
        col_y[n] = str(float(col_y[n]) * pixb)
        col_z[n] = str(float(col_z[n]) * pixc)
    doc.write_file("emda_reboxed_model.cif")    


def get_filename_from_path(file_path):
    import os

    #fpath = os.path.abspath(file_path)
    base=os.path.basename(file_path)
    return os.path.splitext(base)[0]


def _rebox(imap, imask, imodel=None):
    uc, arr, orig = em.get_data(imap)
    _, mask, _ = em.get_data(imask)
    outmapname = get_filename_from_path(imap) + "_emda_reboxed.mrc"
    newarr, newuc = reboxmap(arr, uc, mask=mask, outname=outmapname)
    if imodel is not None:
        outmodelname = get_filename_from_path(imodel) + "_emda_reboxed.cif"
        model_rebox(arr, mask, mmcif_file=imodel, uc=uc)
    return newarr, newuc


def rebox_maps_and_models(maplist, modellist=None, masklist=None):
    from emda.ext.maskmap_class import mask_from_coordinates

    if masklist is not None:
        if modellist is not None:
            assert len(maplist) == len(masklist) == len(modellist)
            # read and rebox
            for imap, imodel, imask in zip(maplist, modellist, masklist):
                _, _ = _rebox(imap=imap, imask=imask, imodel=imodel)
        else:
            assert len(maplist) == len(masklist)
            # read and rebox
            for imap, imask in zip(maplist, masklist):
                _, _ = _rebox(imap=imap, imask=imask)
    else:
        if modellist is not None:
            assert len(modellist) == len(maplist)
            for imap, imodel in zip(maplist, modellist):
                # calculate mask from model each
                _ = mask_from_coordinates(mapname=imap, modelname=imodel)
                # rebox
                _, _ = _rebox(imap=imap, imask='./emda_atomic_mask.mrc')
        else:
            print("Either mask or model is needed for reboxing")
            raise SystemExit("Exiting...")


if __name__=="__main__":
    #imap = "/Users/ranganaw/MRC/REFMAC/EMD-4174/emd_4174.map"
    #imap = "/Users/ranganaw/MRC/REFMAC/ligand_challenge/EMD-7770/Trimmed_maps_and_model/emd_7770_trimmed.mrc"
    #imtz = "/Users/ranganaw/MRC/REFMAC/ligand_challenge/EMD-7770/final_model_map/refined.mtz"
    """ uc, arr, orig = em.get_data(imap)
    map = mtz2map(mtzname=imtz, map_size=arr.shape)
    uc2 = np.array([uc[2], uc[1], uc[0]], 'float')
    em.write_mrc(mapdata=map, filename="refined.mrc", unit_cell=uc2) """
    maplist = [
        "/Users/ranganaw/MRC/REFMAC/COVID19/EMD-11203/emd_11203_half1.map",
        "/Users/ranganaw/MRC/REFMAC/COVID19/EMD-11203/emd_11203_half2.map",
        "/Users/ranganaw/MRC/REFMAC/COVID19/EMD-11203/bestmap1dfsc.mrc"
    ]
    imask = "/Users/ranganaw/MRC/REFMAC/COVID19/EMD-11203/emda_modelmask/emda_atomic_mask.mrc"
    imodel = "/Users/ranganaw/MRC/REFMAC/COVID19/EMD-11203/Refmac5_1/refined.pdb"
    _, mask, _ = em.get_data(imask)
    for imap in maplist:
        uc, arr, orig = em.get_data(imap)
        outname = get_filename_from_path(imap) + "_reboxed.mrc"
        _, _ = reboxmap(arr, uc, mask=mask, outname=outname)

    model_rebox(arr, mask, mmcif_file=imodel, uc=uc)
    



