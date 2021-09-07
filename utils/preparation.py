import numpy as np
import mrcfile
import gemmi
import math


def pad_or_crop_volume(vol, dim_pad=None, pad_value = None, crop_volume=False):
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
            crop_vol = vol[int(round(vol.shape[0]/2-dim_pad[0]/2)):int(round(vol.shape[0]/2+dim_pad[0]/2+dim_pad[0]%2)), :, :]
            crop_vol = crop_vol[:, int(round(vol.shape[1]/2-dim_pad[1]/2)):int(round(vol.shape[1]/2+dim_pad[1]/2+dim_pad[1]%2)), :]
            crop_vol = crop_vol[:, :, int(round(vol.shape[2]/2-dim_pad[2]/2)):int(round(vol.shape[2]/2+dim_pad[2]/2+dim_pad[2]%2))]

            return crop_vol

        else:
            pad_vol = np.pad(vol, ((int(round(dim_pad[0]/2-vol.shape[0]/2)), int(round(dim_pad[0]/2-vol.shape[0]/2+dim_pad[0]%2))), (0,0), (0,0) ), 'constant', constant_values=(pad_value,))
            pad_vol = np.pad(pad_vol, ((0,0), (int(round(dim_pad[1]/2-vol.shape[1]/2)), int(round(dim_pad[1]/2-vol.shape[1]/2+dim_pad[1]%2)) ), (0,0)), 'constant', constant_values=(pad_value,))
            pad_vol = np.pad(pad_vol, ((0,0), (0,0), (int(round(dim_pad[2]/2-vol.shape[2]/2)), int(round(dim_pad[2]/2-vol.shape[2]/2+dim_pad[2]%2)))), 'constant', constant_values=(pad_value,))

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
    #print(mask.shape)
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


def get_spherical_mask(emmap):
    
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
    
def prepare_mask_and_maps_for_scaling(args):
    
    ## Added new part to include pseudo-atomic model
    from pseudomodel_utils import get_modmap_from_pseudomodel
    from emmer.ndimage.map_utils import average_voxel_size, pad_or_crop_image
    from emmer.pdb.pdb_tools import find_wilson_cutoff
    from locscale_headers import run_FDR
    from utils.general import round_up_to_even, round_up_to_odd
    
    emmap = mrcfile.open(args.em_map).data
    verbose = bool(args.verbose)
    
    fsc_resolution = float(args.fsc_resolution)
    
    if args.apix is None:
        apix = average_voxel_size(mrcfile.open(args.em_map).voxel_size)  ## Assuming voxelsize is the same in all directions
    else:
        apix = float(args.apix)
    
    
    ## Get the mask path if provided. Calculate if not provided.
    if args.mask is None:
        if args.verbose:
            print("A mask path has not been provided. False Discovery Rate control (FDR) based confidence map will be calculated at 1% FDR \n")
        if args.fdr_w is None:   # if FDR window size is not set, take window size equal to 10% of emmap height
            fdr_window_size = round_up_to_even(emmap.shape[0] * 0.1)
            print("FDR window size is not set. Using a default window size of {}".format(fdr_window_size))
        else:
            fdr_window_size = int(args.fdr_w)
        
        if args.fdr_filter is not None:
            filter_cutoff = float(args.fdr_filter)
            print("A low pass filter value has been provided. The EM-map will be low pass filtered to {:.2f} A".format(filter_cutoff))
        else:
            filter_cutoff = None
            
        mask_path = run_FDR(emmap_path=args.em_map, window_size = fdr_window_size, fdr=0.01, filter_cutoff=filter_cutoff)
        if mask_path is not None:
            mask = (mrcfile.open(mask_path).data == 1).astype(np.int8)
        else:
            mask = get_spherical_mask(emmap.shape)
    else:
        mask_path = args.mask
        mask = (mrcfile.open(mask_path).data == 1).astype(np.int8)
    
    
    ## Use the mask and emmap to generate a model map using pseudo-atomic model
        
    if args.model_map is None:
        print("Reference Data not supplied! Using pseudo-model")
        pseudomodel_method=args.pm
        pam_distance = float(args.dst)
        if pseudomodel_method is 'random':
            pam_iteration = 100
        elif pseudomodel_method is 'gradient':
            pam_iteration = 50
        
        modmap_path = get_modmap_from_pseudomodel(emmap_path=args.em_map, mask_path=mask_path, 
                                                  pseudomodel_method=pseudomodel_method, pam_distance=pam_distance, pam_iteration=pam_iteration,
                                                  fsc_resolution=fsc_resolution, verbose=verbose)
        modmap = mrcfile.open(modmap_path).data
    else:
        modmap_path = args.model_map
        modmap = mrcfile.open(args.model_map).data        
        
    
    if args.window_size is None:   ## Use default window size of 25 A
        wn = round_up_to_odd(25 / apix)+1
        if verbose:
            print("Using a default window size of 25 A, corresponding to approximately {} pixels".format(wn))
        
    elif args.window_size is not None:
        wn = round_up_to_odd(float(args.window_size))
        if verbose:
            print("Provided window size in pixels is {} corresponding to {:.2f} Angstorm".format(wn, wn*apix))

    window_bleed_and_pad = check_for_window_bleeding(mask, wn)
    if window_bleed_and_pad:
        pad_int_emmap = compute_padding_average(emmap, mask)
        pad_int_modmap = compute_padding_average(modmap, mask)
        map_shape = [(emmap.shape[0] + wn), (emmap.shape[1] + wn), (emmap.shape[2] + wn)]
        emmap = pad_or_crop_volume(emmap, map_shape, pad_int_emmap)
        modmap = pad_or_crop_volume(modmap, map_shape, pad_int_modmap)
        mask = pad_or_crop_volume(mask, map_shape, 0)
    
    ## 
    wilson_cutoff = find_wilson_cutoff(mask_path=mask_path)
    fsc_cutoff = fsc_resolution
    if verbose:
        print("Using Wilson Cutoff of: {} and FSC cutoff of {}".format(wilson_cutoff, fsc_cutoff))
    frequency_cutoffs = (wilson_cutoff, fsc_cutoff)
    
    return emmap, modmap, mask, wn, window_bleed_and_pad, apix, frequency_cutoffs
