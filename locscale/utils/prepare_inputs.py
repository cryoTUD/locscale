import numpy as np
import mrcfile
import locscale.include.emmer as emmer


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

def round_up_proper(x):
    epsilon = 1e-5  ## To round up in case of rounding to odd
    return np.round(x+epsilon).astype(int)

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


def generate_filename_from_halfmap_path(in_path):
    ## find filename in the path    
    filename = in_path.split("/")[-1]
    
    ## find EMDB ID in filename
    
    possible_emdb_id = [filename[x:x+4] for x in range(len(filename)-3) if filename[x:x+4].isnumeric()]
    if len(possible_emdb_id) == 1:
        emdb_id = possible_emdb_id[0]
        newfilename = ["EMD_"+emdb_id+"_unfiltered.mrc"]
    else:
        newfilename = ["emdb_map_unfiltered.mrc"]
    
    
    new_path = "/".join(in_path.split("/")[:-1]+newfilename)
    
    return new_path
    
    
def prepare_mask_and_maps_for_scaling(args):
    '''
    Parse the command line arguments and return inputs for computing local amplitude scaling 

    Parameters
    ----------
    args : Namespace

    Returns
    -------
    parsed_inputs_dict : dict
        Parsed inputs dictionary

    '''

    print("Preparing your inputs for LocScale")
    from locscale.pseudomodel.pipeline import get_modmap
    from locscale.pseudomodel.pseudomodel_headers import number_of_segments, run_FDR, run_mapmask, check_dependencies
    from locscale.utils.general import round_up_to_even, round_up_to_odd, get_emmap_path_from_args
    from locscale.include.emmer.ndimage.map_tools import add_half_maps
    from locscale.include.emmer.ndimage.map_utils import average_voxel_size
    from locscale.include.emmer.ndimage.profile_tools import compute_radial_profile, estimate_bfactor_through_pwlf, frequency_array
    from locscale.include.emmer.pdb.pdb_tools import find_wilson_cutoff
    from locscale.include.emmer.pdb.pdb_utils import shift_coordinates
    
    print("Check relevant paths for LocScale \n")
    print(check_dependencies())
    
    scale_using_theoretical_profile = not(args.ignore_profiles)
    
    emmap_path, shift_vector = get_emmap_path_from_args(args)
    
    xyz_emmap_path = run_mapmask(emmap_path)  
    
    
    ## run_mapmask() function makes the axis order of the .mrc file to XYZ 
    
    xyz_emmap = mrcfile.open(xyz_emmap_path).data
    
    verbose = bool(args.verbose)
    
    fsc_resolution = float(args.ref_resolution)
    
    if args.apix is None:
        apix = average_voxel_size(mrcfile.open(emmap_path).voxel_size)  ## Assuming voxelsize is the same in all directions
    else:
        apix = float(args.apix)
    
    
    ## Get the mask path if provided. Calculate if not provided.
    
    if args.mask is None:
        if args.verbose:
            print("A mask path has not been provided. False Discovery Rate control (FDR) based confidence map will be calculated at 1% FDR \n")
        if args.fdr_window_size is None:   # if FDR window size is not set, take window size equal to 10% of emmap height
            fdr_window_size = round_up_to_even(xyz_emmap.shape[0] * 0.1)
            print("FDR window size is not set. Using a default window size of {} \n".format(fdr_window_size))
        else:
            fdr_window_size = int(args.fdr_w)
        
        if args.fdr_filter is not None:
            filter_cutoff = float(args.fdr_filter)
            print("A low pass filter value has been provided. The EM-map will be low pass filtered to {:.2f} A \n".format(filter_cutoff))
        else:
            filter_cutoff = None
            
        mask_path = run_FDR(emmap_path=emmap_path, window_size = fdr_window_size, fdr=0.01, filter_cutoff=filter_cutoff)
        xyz_mask_path = run_mapmask(mask_path)
        
        if xyz_mask_path is not None:
            xyz_mask = (mrcfile.open(xyz_mask_path).data == 1).astype(np.int8)
        else:
            xyz_mask = get_spherical_mask(xyz_emmap.shape)
    else:
        mask_path = args.mask
        xyz_mask_path = run_mapmask(mask_path)
        xyz_mask = (mrcfile.open(xyz_mask_path).data == 1).astype(np.int8)
    
    
    ## Use the mask and emmap to generate a model map using pseudo-atomic model
    if args.molecular_weight is not None:
        molecular_weight = float(args.molecular_weight)    
    else:
        molecular_weight = None

    if args.model_map is None:
        
        pdb_path = args.model_coordinates
        if pdb_path is not None:
            scale_using_theoretical_profile = False ## If a PDB_path is provided, assume that it is an atomic model thus set this flag as False
            shift_coordinates(in_model_path=pdb_path, trans_matrix=shift_vector,
                                         out_model_path=pdb_path[:-4]+"_shifted.pdb")
            pdb_path = pdb_path[:-4]+"_shifted.pdb"
            
        add_blur = float(args.add_blur)
        skip_refine = args.skip_refine
        model_resolution = args.model_resolution
        ## Defaults for pseudo-atomic model 
        pseudomodel_method=args.pseudomodel_method
        pam_distance = float(args.distance)
        refmac_iter = int(args.refmac_iterations)
        if pseudomodel_method == 'random' and args.total_iterations is None:
            pam_iteration = 100
        elif pseudomodel_method == 'gradient' and args.total_iterations is None:
            pam_iteration = 50
        elif args.total_iterations is not None:
            pam_iteration = int(args.total_iterations)
        build_ca_only = args.build_ca_only
        
        ## Get reference map using get_modmap_from_pseudomodel()
        ## Note that if a pdb_path is provided then the function 
        ## will use that instead of running pseudo-atomic model 
        ## routine. 
        
        modmap_args = {
            'emmap_path':xyz_emmap_path,
            'mask_path':xyz_mask_path,
            'pdb_path':pdb_path,
            'pseudomodel_method':pseudomodel_method,
            'pam_distance':pam_distance,
            'pam_iteration':pam_iteration,
            'fsc_resolution':fsc_resolution,
            'refmac_iter':refmac_iter,
            'add_blur':add_blur,
            'skip_refine':skip_refine,
            'model_resolution':model_resolution,
            'pg_symmetry':args.symmetry,
            'molecular_weight':molecular_weight,
            'build_ca_only':build_ca_only,
            'verbose':verbose
        }
        
        modmap_path = get_modmap(modmap_args)
        
        xyz_modmap_path = run_mapmask(modmap_path, return_same_path=True)
        xyz_modmap = mrcfile.open(xyz_modmap_path).data
    else:
        scale_using_theoretical_profile = False ## If a model map is provided, assume that it is from an atomic model thus set this flag as False no matter what the user input 
        modmap_path = args.model_map
        xyz_modmap_path = run_mapmask(modmap_path)
        xyz_modmap = mrcfile.open(xyz_modmap_path).data        
        
    
    if args.window_size is None:   ## Use default window size of 25 A
        wn = round_up_to_even(25 / apix)
        if verbose:
            print("Using a default window size of 25 A, corresponding to approximately {} pixels".format(wn))
        
    elif args.window_size is not None:
        wn = round_up_to_even(int(args.window_size))
        if verbose:
            print("Provided window size in pixels is {} corresponding to {:.2f} Angstorm".format(wn, wn*apix))

    window_bleed_and_pad = check_for_window_bleeding(xyz_mask, wn)
    if window_bleed_and_pad:
        pad_int_emmap = compute_padding_average(xyz_emmap, xyz_mask)
        pad_int_modmap = compute_padding_average(xyz_modmap, xyz_mask)
        map_shape = [(xyz_emmap.shape[0] + wn), (xyz_emmap.shape[1] + wn), (xyz_emmap.shape[2] + wn)]
        xyz_emmap = pad_or_crop_volume(xyz_emmap, map_shape, pad_int_emmap)
        xyz_modmap = pad_or_crop_volume(xyz_modmap, map_shape, pad_int_modmap)
        xyz_mask = pad_or_crop_volume(xyz_mask, map_shape, 0)
    

    ## Next few lines of code characterizes radial profile of 
    ## input emmap : 
        ## wilson cutoff : threshold between guinier and wilson regimes in the radial profile
        ## high frequency cutoff : threshold above which to computing bfactor becomes valid  (for low resolution map, it's same as wilson cutoff)
        ## FSC cutoff : threshold above which amplitudes of signal becomes weaked compared to noise
        
    
    wilson_cutoff = find_wilson_cutoff(mask_path=xyz_mask_path)
    smooth_factor = args.smooth_factor
    boost_secondary_structure = args.boost_secondary_structure
    if fsc_resolution > 6:
        high_frequency_cutoff = wilson_cutoff
        nyquist = (round(2*apix*10)+1)/10
        #fsc_cutoff = fsc_resolution
        bfactor_info = [0,np.array([0,0,0]),np.array([0,0,0])]
        pwlf_fit_quality = 0
    else:
        rp_emmap = compute_radial_profile(xyz_emmap)
        freq = frequency_array(amplitudes=rp_emmap, apix=apix)
        num_segments = number_of_segments(fsc_resolution)
        bfactor, amp, (fit,z,slope) = estimate_bfactor_through_pwlf(freq=freq, amplitudes=rp_emmap, wilson_cutoff=wilson_cutoff, fsc_cutoff=fsc_resolution,num_segments=num_segments)
        nyquist = (round(2*apix*10)+1)/10
        #fsc_cutoff = fsc_resolution
        high_frequency_cutoff = 1/np.sqrt(z[-2])
        bfactor_info = [round(bfactor,2), 1/np.sqrt(z).round(2), np.array(slope).round(2)]  ## For information at end
        pwlf_fit_quality = fit.r_squared()
    
    dev_mode = args.dev_mode
    if dev_mode:
        print("DEV MODE: Scaling using theoretical profiles even if you have input an atomic model!")
        print("Before: scale_using_theoretical_profile=",scale_using_theoretical_profile)
        scale_using_theoretical_profile = None
        scale_using_theoretical_profile = True
        print("After: scale_using_theoretical_profile=",scale_using_theoretical_profile)        
    
    if verbose and scale_using_theoretical_profile:
        print("To compute bfactors of local windows: \nUsing High Frequency Cutoff of: {:.2f} and FSC cutoff of {}".format(high_frequency_cutoff, nyquist))
        print("To merge reference and theoretical profiles: \n")
        print("Using Wilson cutoff of {:.2f} A and smooth factor of {:.2f}".format(wilson_cutoff, smooth_factor))
        
    
    
    scale_factor_arguments = {}
    scale_factor_arguments['wilson'] = wilson_cutoff
    scale_factor_arguments['high_freq'] = high_frequency_cutoff
    scale_factor_arguments['fsc_cutoff'] = fsc_resolution
    scale_factor_arguments['nyquist'] = nyquist
    scale_factor_arguments['smooth'] = smooth_factor
    scale_factor_arguments['boost_secondary_structure'] = boost_secondary_structure
    
    if verbose:
        
        print("Preparation completed. Now running LocScale!")
    
    
    parsed_inputs_dict = {}
    parsed_inputs_dict['emmap'] = xyz_emmap
    parsed_inputs_dict['modmap'] = xyz_modmap
    parsed_inputs_dict['mask'] = xyz_mask
    parsed_inputs_dict['wn'] = wn
    parsed_inputs_dict['apix'] = apix
    parsed_inputs_dict['use_theoretical'] = scale_using_theoretical_profile
    parsed_inputs_dict['scale_factor_args'] = scale_factor_arguments
    parsed_inputs_dict['verbose'] = verbose
    parsed_inputs_dict['win_bleed_pad'] = window_bleed_and_pad
    parsed_inputs_dict['bfactor_info'] = bfactor_info
    parsed_inputs_dict['fsc_resolution'] = fsc_resolution
    parsed_inputs_dict['PWLF_fit'] = pwlf_fit_quality
    parsed_inputs_dict['emmap_path'] = xyz_emmap_path
    parsed_inputs_dict['mask_path'] = xyz_mask_path
    
    ## all maps should have same shape
    assert xyz_emmap.shape == xyz_modmap.shape == xyz_mask.shape, "The input maps and mask do not have the same shape"
    ## emmap and modmap should not be zeros
    assert abs(xyz_emmap.sum()) > 0 and abs(xyz_modmap.sum()) > 0, "Emmap and Modmap should not be zeros!"
    ## No element of the mask should be negative
    assert (xyz_mask>=0).any(), "Negative numbers found in mask"
    
    
    
    
    
    
    return parsed_inputs_dict
