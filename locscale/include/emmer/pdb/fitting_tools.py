# Set of tools related to the fitting of an atomic model to a density map
import numpy as np
from locscale.include.emmer.ndimage.map_utils import load_map
from locscale.include.emmer.pdb.pdb_to_map import detect_pdb_input

def map_values_at_atomic_locations(input_pdb, emmap_path, sample_percentage=None):
    '''
    Returns a dictionary which contains the value of the map at the location of each atom in the pdb file
    map_values = dict
    key = atomic coordinate (XYZ, in Angstroms)
    value = map value at that coordinate

    '''
    from locscale.include.emmer.pdb.pdb_utils import get_coordinates, get_atomic_point_map
    from locscale.include.emmer.ndimage.map_utils import convert_mrc_to_pdb_position, convert_pdb_to_mrc_position

    # Load the inputs
    st = detect_pdb_input(input_pdb)
    emmap, apix = load_map(emmap_path)

    # Get the coordinates of the atoms in the pdb file
    coords = get_coordinates(st)

    # Convert the coordinates to the map coordinates
    mrc_coords = np.array(convert_pdb_to_mrc_position(coords, apix))

    atomic_point_map = get_atomic_point_map(mrc_coords, emmap.shape).astype(bool)
    map_values = emmap[atomic_point_map]
    return map_values
    

    # # sample the indices if needed
    # if sample_percentage is not None:
    #     assert sample_percentage <= 1 and sample_percentage > 0, 'sample_percentage must be between 0 and 1'
    #     n = len(mrc_coords)
    #     sample_size = int(n * sample_percentage)
    #     random_indices = np.random.choice(n, sample_size, replace=False).astype(int)
    #     mrc_coords_sample = mrc_coords[random_indices]
    # else:
    #     mrc_coords_sample = mrc_coords


    # # Get the values of the map at the coordinates of the atoms
    # map_values = {}
    # for i, mrc_coord in enumerate(mrc_coords_sample):
    #     map_values[tuple(coords[i])] = emmap[mrc_coord[0], mrc_coord[1], mrc_coord[2]]
    

    # return map_values


def translation_test_pdb(pdb_path, emmap_path, sample_percentage=None):
    from locscale.include.emmer.ndimage.map_utils import load_map
    from locscale.include.emmer.pdb.pdb_to_map import detect_pdb_input
    from locscale.include.emmer.pdb.pdb_utils import shift_coordinates
    from tqdm import tqdm

    # Load the inputs
    st = detect_pdb_input(pdb_path)

    # define translations array
    translations_magnitudes = np.linspace(-5,5,10)

    # apply shifts to the pdb
    shifted_mean_map_values_x = {}
    shifted_mean_map_values_y = {}
    shifted_mean_map_values_z = {}

    for d in tqdm(translations_magnitudes, desc='Translation test'):
        translation_matrix_x = [0,0,d]
        translation_matrix_y = [0,d,0]
        translation_matrix_z = [d,0,0]

        shifted_structure_x = shift_coordinates(input_structure=st, trans_matrix=translation_matrix_x)
        shifted_structure_y = shift_coordinates(input_structure=st, trans_matrix=translation_matrix_y)
        shifted_structure_z = shift_coordinates(input_structure=st, trans_matrix=translation_matrix_z)

        shifted_map_values_x = map_values_at_atomic_locations(shifted_structure_x, emmap_path, sample_percentage=sample_percentage)
        shifted_map_values_y = map_values_at_atomic_locations(shifted_structure_y, emmap_path, sample_percentage=sample_percentage)
        shifted_map_values_z = map_values_at_atomic_locations(shifted_structure_z, emmap_path, sample_percentage=sample_percentage)

        mean_shifted_map_values_x = shifted_map_values_x.mean()
        mean_shifted_map_values_y = shifted_map_values_y.mean()
        mean_shifted_map_values_z = shifted_map_values_z.mean()

        shifted_mean_map_values_x[d] = mean_shifted_map_values_x
        shifted_mean_map_values_y[d] = mean_shifted_map_values_y
        shifted_mean_map_values_z[d] = mean_shifted_map_values_z

    shifted_mean_values = {
        'x': shifted_mean_map_values_x,
        'y': shifted_mean_map_values_y,
        'z': shifted_mean_map_values_z
    }

    return shifted_mean_values

def check_pdb_fit_with_map(pdb_path, emmap_path):

    shifted_mean_values = translation_test_pdb(pdb_path, emmap_path)

    translation_magnitudes = np.array(list(shifted_mean_values['x'].keys()))
    mean_shifted_map_values_x = np.array(list(shifted_mean_values['x'].values()))
    mean_shifted_map_values_y = np.array(list(shifted_mean_values['y'].values()))
    mean_shifted_map_values_z = np.array(list(shifted_mean_values['z'].values()))

    

def compute_model_to_map_correlation(pdb_path, emmap_path):
    from locscale.include.emmer.ndimage.map_utils import load_map
    from locscale.include.emmer.ndimage.map_tools import get_atomic_model_mask, compute_real_space_correlation
    from locscale.include.emmer.pdb.pdb_to_map import detect_pdb_input, pdb2map

    atomic_model_mask = get_atomic_model_mask(emmap_path=emmap_path, pdb_path=pdb_path, save_files=False)
    atomic_mask_binarised = (atomic_model_mask > 0.5).astype(bool)
    emmap, apix = load_map(emmap_path)

    simmap = pdb2map(pdb_path, apix=apix, size=emmap.shape)

    # compute the correlation
    rscc = compute_real_space_correlation(emmap[atomic_mask_binarised], simmap[atomic_mask_binarised])

    return rscc







