import numpy as np
import math

## MAP PROCESSING FUNCTIONS REQUIRED FOR LOCSCALE 

def shift_map_to_zero_origin(emmap_path):
    '''
    Determines the map origin from header file and changes it to zero

    Parameters
    ----------
    emmap_path : str
        DESCRIPTION.

    Returns
    -------
    shift_vector : numpy.ndarray (len=3)

    '''    
    import mrcfile
    from locscale.include.emmer.ndimage.map_utils import save_as_mrc
    
    target_origin = np.array([0,0,0])
    voxel_size = np.array(mrcfile.open(emmap_path).voxel_size.tolist())
    current_origin = np.array(mrcfile.open(emmap_path).header.origin.tolist()) 
    
    emmap_data = mrcfile.open(emmap_path).data
    
    output_file = emmap_path
    save_as_mrc(map_data=emmap_data, output_filename=emmap_path, apix=voxel_size, origin=0)
    
    shift_vector = target_origin - current_origin
    return shift_vector

def get_spherical_mask(emmap):
    from locscale.utils.general import pad_or_crop_volume
    
    mask = np.zeros(emmap.shape)

    if mask.shape[0] == mask.shape[1] and mask.shape[0] == mask.shape[2] and mask.shape[1] == mask.shape[2]:
        rad = mask.shape[0] // 2
        z,y,x = np.ogrid[-rad: rad+1, -rad: rad+1, -rad: rad+1]
        mask = (x**2+y**2+z**2 <= rad**2).astype(np.int_).astype(np.int8)
        mask = pad_or_crop_volume(mask,emmap.shape)
        mask = (mask==1).astype(np.int8)
    else:
        mask += 1
        mask = mask[0:mask.shape[0]-1, 0:mask.shape[1]-1, 0:mask.shape[2]-1]
        mask = pad_or_crop_volume(emmap, (emmap.shape), pad_value=0)
    
    return mask

def put_scaled_voxels_back_in_original_volume_including_padding(sharpened_vals, masked_indices, map_shape):
    map_scaled = np.zeros(np.prod(map_shape))
    map_scaled[masked_indices] = sharpened_vals
    map_scaled = map_scaled.reshape(map_shape)

    return map_scaled

def pad_or_crop_volume(vol, dim_pad=None, pad_value = None, crop_volume=False):
    from locscale.utils.math_tools import round_up_proper
    if (dim_pad == None):
        return vol
    else:
        dim_pad = np.round(np.array(dim_pad)).astype('int')
        #print(dim_pad)

        if pad_value == None:
            pad_value = 0

        if (dim_pad[0] <= vol.shape[0] or dim_pad[1] <= vol.shape[1] or dim_pad[2] <= vol.shape[2]):
            crop_volume = True

        if crop_volume:
            k_start = round_up_proper(vol.shape[0]/2-dim_pad[0]/2)
            k_end = round_up_proper(vol.shape[0]/2+dim_pad[0]/2)
            j_start = round_up_proper(vol.shape[1]/2-dim_pad[1]/2)
            j_end = round_up_proper(vol.shape[1]/2+dim_pad[1]/2)
            i_start = round_up_proper(vol.shape[2]/2-dim_pad[2]/2)
            i_end = round_up_proper(vol.shape[2]/2+dim_pad[2]/2)
            crop_vol = vol[k_start:k_end, :, :]
            crop_vol = crop_vol[:, j_start:j_end, :]
            crop_vol = crop_vol[:, :, i_start:i_end]

            return crop_vol

        else:
            k_start = round_up_proper(dim_pad[0]/2-vol.shape[0]/2)
            k_end = round_up_proper(dim_pad[0]/2-vol.shape[0]/2)
            j_start = round_up_proper(dim_pad[1]/2-vol.shape[1]/2)
            j_end = round_up_proper(dim_pad[1]/2-vol.shape[1]/2)
            i_start = round_up_proper(dim_pad[2]/2-vol.shape[2]/2)
            i_end = round_up_proper(dim_pad[2]/2-vol.shape[2]/2)
            
            pad_vol = np.pad(vol, ((k_start, k_end ), (0,0), (0,0) ), 'constant', constant_values=(pad_value,))
            pad_vol = np.pad(pad_vol, ((0,0), (j_start, j_end ), (0,0)), 'constant', constant_values=(pad_value,))
            pad_vol = np.pad(pad_vol, ((0,0), (0,0), (i_start, i_end )), 'constant', constant_values=(pad_value,))
            
            if pad_vol.shape[0] != dim_pad[0] or pad_vol.shape[1] != dim_pad[1] or pad_vol.shape[2] != dim_pad[2]:
                print("Requested pad volume shape {} not equal to the shape of the padded volume returned{}. Input map shape might be an odd sized map.".format(dim_pad, pad_vol.shape))
                

            return pad_vol

def compute_padding_average(vol, mask):
    mask = (mask == 1).astype(np.int8)
    #inverted_mask = np.logical_not(mask)
    average_padding_intensity = np.mean(np.ma.masked_array(vol, mask))
    return average_padding_intensity

def get_xyz_locs_and_indices_after_edge_cropping_and_masking(mask, wn):
    mask = np.copy(mask)
    nk, nj, ni = mask.shape

    kk, jj, ii = np.indices((mask.shape))
    kk_flat = kk.ravel()
    jj_flat = jj.ravel()
    ii_flat = ii.ravel()

    mask_bin = np.array(mask.ravel(), dtype=np.bool)
    indices = np.arange(mask.size)
    masked_indices = indices[mask_bin]
    cropped_indices = indices[(wn / 2 <= kk_flat) & (kk_flat < (nk - wn / 2)) &
                              (wn / 2 <= jj_flat) & (jj_flat < (nj - wn / 2)) &
                              (wn / 2 <= ii_flat) & (ii_flat < (ni - wn / 2))]

    cropp_n_mask_ind = np.intersect1d(masked_indices, cropped_indices)

    xyz_locs = np.column_stack((kk_flat[cropp_n_mask_ind], jj_flat[cropp_n_mask_ind], ii_flat[cropp_n_mask_ind]))

    return xyz_locs, cropp_n_mask_ind, mask.shape

def check_for_window_bleeding(mask,wn):
    from locscale.utils.general import get_xyz_locs_and_indices_after_edge_cropping_and_masking
    masked_xyz_locs, masked_indices, mask_shape = get_xyz_locs_and_indices_after_edge_cropping_and_masking(mask, 0)

    zs, ys, xs = masked_xyz_locs.T
    nk, nj, ni = mask_shape
    #print(xs.shape, ys.shape, zs.shape)
    #print(nk,nj,ni)
    #print(wn)

    if xs.min() < wn / 2 or xs.max() > (ni - wn / 2) or \
    ys.min() < wn / 2 or ys.max() > (nj - wn / 2) or \
    zs.min() < wn / 2 or zs.max() > (nk - wn / 2):
        window_bleed = True
    else:
        window_bleed = False

    return window_bleed

def normalise_intensity_levels(from_emmap, to_levels=[0,1]):
    normalise_between_zero_one = (from_emmap - from_emmap.min()) / (from_emmap.max() - from_emmap.min())
    to_levels = np.array(to_levels)
    
    min_value = to_levels.min()
    max_value = to_levels.max()
    scale_factor = max_value-min_value
    
    normalised = min_value + normalise_between_zero_one * scale_factor
    
    return normalised    

def shift_radial_profile(from_emmap, to_emmap):
    '''
    To shift the radial profile of one emmap so that DC power matches another emmap

    Parameters
    ----------
    from_emmap : numpy.ndarray
        DESCRIPTION.
    to_emmap : numpy.ndarray
        DESCRIPTION.

    Returns
    -------
    emmap_shifted_rp : numpy.ndarray

    '''
    from locscale.include.emmer.ndimage.profile_tools import compute_radial_profile, plot_radial_profile, frequency_array
    from locscale.include.emmer.ndimage.map_tools import set_radial_profile, compute_scale_factors
    
    rp_from_emmap = compute_radial_profile(from_emmap)
    rp_to_emmap, radii = compute_radial_profile(to_emmap, return_indices=True)
    
    DC_power_diff = rp_to_emmap.max() - rp_from_emmap.max()
    print((DC_power_diff))
    new_rp_from_emmap = rp_from_emmap * 20
    scale_factors = compute_scale_factors(rp_from_emmap, new_rp_from_emmap)
    emmap_shifted_rp = set_radial_profile(from_emmap, scale_factors, radii)
    rp_after_shifted = compute_radial_profile(emmap_shifted_rp)
    freq = frequency_array(rp_from_emmap, 1.2156)
    plot_radial_profile(freq,[rp_from_emmap, rp_to_emmap, new_rp_from_emmap,rp_after_shifted])
    
    
    return emmap_shifted_rp
    
##### SAVE FUNCTIONS #####

def save_list_as_map(values_list, masked_indices, map_shape, map_path, apix):
    from locscale.include.emmer.ndimage.map_utils import save_as_mrc
    from locscale.utils.general import put_scaled_voxels_back_in_original_volume_including_padding
    value_map = put_scaled_voxels_back_in_original_volume_including_padding(values_list, masked_indices, map_shape)
    save_as_mrc(value_map, output_filename=map_path, apix=apix)

def write_out_final_volume_window_back_if_required(args, LocScaleVol, parsed_inputs_dict):
    from locscale.utils.general import pad_or_crop_volume
    from locscale.include.emmer.ndimage.map_utils import save_as_mrc
    from locscale.utils.plot_tools import make_locscale_report
    import mrcfile
    
    input_map = mrcfile.open(parsed_inputs_dict['emmap_path']).data
    
    wn = parsed_inputs_dict['wn']
    window_bleed_and_pad =parsed_inputs_dict['win_bleed_pad']
    apix = parsed_inputs_dict['apix']
        
    if window_bleed_and_pad:
        #map_shape = [(LocScaleVol.shape[0] - wn), (LocScaleVol.shape[1] - wn), (LocScaleVol.shape[2] - wn)]
        map_shape = input_map.shape
        LocScaleVol = pad_or_crop_volume(LocScaleVol, (map_shape))
    output_filename = args.outfile
    if args.dev_mode:
        output_filename = output_filename[:-4]+"_devmode.mrc"
    save_as_mrc(map_data=LocScaleVol, output_filename=output_filename, apix=apix, origin=0, verbose=True)
    
    if args.symmetry != "C1":
        resolution = parsed_inputs_dict['fsc_resolution']
        print("Imposing a symmetry condition of {}".format(args.symmetry))
        from locscale.include.symmetry_emda.symmetry_map import symmetrize_map_emda
        
        LocScaleVol_sym = symmetrize_map_emda(emmap_path=output_filename,pg=args.symmetry)

        output_filename = output_filename[:-4]+"_symmetrised.mrc"

        save_as_mrc(map_data=LocScaleVol_sym, output_filename=output_filename, apix=apix, origin=0, verbose=True)
    
    if args.output_report:
        make_locscale_report(args, parsed_inputs_dict, output_filename, window_bleed_and_pad)
    
    return LocScaleVol

##### MPI related functions #####

def split_sequence_evenly(seq, size):
    """
    >>> split_sequence_evenly(list(range(9)), 4)
    [[0, 1], [2, 3, 4], [5, 6], [7, 8]]
    >>> split_sequence_evenly(list(range(18)), 4)
    [[0, 1, 2, 3, 4], [5, 6, 7, 8], [9, 10, 11, 12, 13], [14, 15, 16, 17]]
    """
    from locscale.utils.math_tools import round_up_proper
    newseq = []
    splitsize = 1.0 / size * len(seq)
    for i in range(size):
        newseq.append(seq[round_up_proper(i * splitsize):round_up_proper((i+1) * splitsize)])
    return newseq

def merge_sequence_of_sequences(seq):
    """
    >>> merge_sequence_of_sequences([list(range(9)), list(range(3))])
    [0, 1, 2, 3, 4, 5, 6, 7, 8, 0, 1, 2]
    >>> merge_sequence_of_sequences([list(range(9)), [], list(range(3))])
    [0, 1, 2, 3, 4, 5, 6, 7, 8, 0, 1, 2]
    """
    newseq = [number for sequence in seq for number in sequence]

    return newseq
