import gemmi
import numpy as np

def center_of_mass_model(mmcif_file):
    st = gemmi.read_structure(mmcif_file)
    model = st[0]
    com = model.calculate_center_of_mass()
    return np.array([com.x, com.y, com.z], 'float')

def shift_to_origin(mmcif_file):
    doc = gemmi.cif.read_file(mmcif_file)
    st = gemmi.read_structure(mmcif_file)
    model = st[0]
    com = model.calculate_center_of_mass()
    block = doc.sole_block()
    col_x = block.find_values("_atom_site.Cartn_x")
    col_y = block.find_values("_atom_site.Cartn_y")
    col_z = block.find_values("_atom_site.Cartn_z")
    vec = np.zeros(3, dtype='float')
    for n, _ in enumerate(col_x):
        col_x[n] = str(float(col_x[n]) - com.x)
        col_y[n] = str(float(col_y[n]) - com.y)
        col_z[n] = str(float(col_z[n]) - com.z)
    return doc

def shift_model(mmcif_file, shift):
    doc = gemmi.cif.read_file(mmcif_file)
    block = doc.sole_block()
    col_x = block.find_values("_atom_site.Cartn_x")
    col_y = block.find_values("_atom_site.Cartn_y")
    col_z = block.find_values("_atom_site.Cartn_z")
    vec = np.zeros(3, dtype='float')
    for n, _ in enumerate(col_x):
        col_x[n] = str(float(col_x[n]) + shift[0])
        col_y[n] = str(float(col_y[n]) + shift[1])
        col_z[n] = str(float(col_z[n]) + shift[2])
    return doc