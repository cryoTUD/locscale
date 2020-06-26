import mrcfile

import numpy as np

from Bio.PDB.PDBParser import PDBParser

from Bio.PDB.MMCIFParser import MMCIFParser

from Bio.PDB.PDBIO import PDBIO

from Bio.PDB.mmcifio import MMCIFIO

import argparse





def prepare_locscale_map_model(in_model, in_map, out_model, out_map):

    """

    Sets nxstart, nystart and nzstart to 0, applies the same translation to the atomic model.

    Applies padding to the grid in order to avoid grid size prime factors greater than 19.

    """

    with mrcfile.mmap(in_map, mode='r+', permissive=True) as mrc:

        voxel_size = (mrc.header.cella.x / mrc.header.nx,

                      mrc.header.cella.y / mrc.header.ny,

                      mrc.header.cella.z / mrc.header.nz)

        xyz_start = (mrc.header.nxstart, mrc.header.nystart, mrc.header.nzstart)

        xyz = mrc.data

        cell = (mrc.header.cella.x, mrc.header.cella.y, mrc.header.cella.z)



    grid = xyz.shape

    new_grid = find_good_grid(grid)

    if new_grid != grid:

        pad = []

        for i in range(3):

            pad.append(new_grid[i] - grid[i])



        new_cell = (cell[0] + voxel_size[0] * pad[0],

                    cell[1] + voxel_size[1] * pad[1],

                    cell[2] + voxel_size[2] * pad[2])

        cell = new_cell

        xyz = np.pad(xyz, ((0, pad[0]), (0, pad[1]), (0, pad[2])), mode='constant', constant_values=0)



    with mrcfile.new(out_map, overwrite=True) as mr:

        mr.set_data(xyz)

        mr.update_header_from_data()

        mr.header.cella = cell

        mr.header.nxstart = 0

        mr.header.nystart = 0

        mr.header.nzstart = 0



    if xyz_start != (0, 0, 0):

        tr_matrix = []

        for i in (2, 1, 0):

            tr_matrix.append(xyz_start[i]*voxel_size[i])

    else:

        tr_matrix = [0, 0, 0]



    shift_coordinates(in_model, out_model, tr_matrix)





def largest_prime_factor(n):

    i = 2

    while i * i <= n:

        if n % i:

            i += 1

        else:

            n //= i

    return n





def find_good_grid(in_grid):

    out_grid = []

    for dimension in in_grid:

        new_dimension = dimension

        if new_dimension % 2:

            new_dimension += 1

        largest_prime = largest_prime_factor(new_dimension)

        while largest_prime > 19:

            new_dimension += 2

            largest_prime = largest_prime_factor(new_dimension)

        out_grid.append(new_dimension)

    return tuple(out_grid)





def shift_coordinates(in_model_path, out_model_path, trans_matrix):

    if in_model_path.split('.')[-1] == 'pdb':

        file_type = 'pdb'

        parser = PDBParser(PERMISSIVE=1)

        structure = parser.get_structure('in_model', in_model_path)

    elif in_model_path.split('.')[-1] == 'cif':

        file_type = 'cif'

        parser = MMCIFParser()

        structure = parser.get_structure('in_model', in_model_path)

    else:

        raise Exception('Please provide a the input model in PDB or CIF format')

    for atom in structure.get_atoms():

        atom.get_coord()[0] = atom.get_coord()[0] - trans_matrix[0]

        atom.get_coord()[1] = atom.get_coord()[1] - trans_matrix[1]

        atom.get_coord()[2] = atom.get_coord()[2] - trans_matrix[2]



    if file_type == 'pdb:':

        io = PDBIO()

    else:

        io = MMCIFIO()

    io.set_structure(structure)

    io.save(out_model_path)





def main():

    parser = argparse.ArgumentParser(description='Prepare the map and atomic model for ccpem-locscale')

    parser.add_argument("-i_mo", "--in_model", dest="in_model", required=True, help="Input atomic model file")

    parser.add_argument("-i_ma", "--in_map", dest="in_map", required=True, help="Input density map file")

    parser.add_argument("-o_mo", "--out_model", dest="out_model", required=True, help="Output atomic model file")

    parser.add_argument("-o_ma", "--out_map", dest="out_map", required=True, help="Output density map file")

    args = parser.parse_args()



    prepare_locscale_map_model(args.in_model, args.in_map, args.out_model, args.out_map)





if __name__ == '__main__':

    main()
