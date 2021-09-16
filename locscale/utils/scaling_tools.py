import numpy as np
#import gemmi


def get_theoretical_profile(length,apix):
    import pickle
    from locscale.include.emmer.ndimage.profile_tools import resample_1d
    from locscale.pseudomodel.pseudomodel_headers import check_dependencies
    
    
    path_to_locscale = check_dependencies()['locscale']
    location_of_theoretical_profiles = path_to_locscale + "/locscale/utils/theoretical_profiles.pickle"
    
    with open(location_of_theoretical_profiles,'rb') as f:
        profiles = pickle.load(f)
    
    frequency_limits = (float(1/(apix*length)),float(1/(apix*2)))
    helix_profile = profiles['helix']
    resampled_helix_profile = resample_1d(helix_profile['freq'], helix_profile['profile'],num=length,xlims=frequency_limits)
    return resampled_helix_profile


def compute_radial_profile(vol, center=[0,0,0], return_indices=False):
    dim = vol.shape
    m = np.mod(vol.shape,2)
    # make compliant with both fftn and rfftn
    if center is None:
        ps = np.abs(np.fft.fftshift((np.fft.fftn(vol))))
        z, y, x = np.indices(ps.shape)
        center = tuple((a - 1) / 2.0 for a in ps.shape[::-1])
        radii = np.sqrt((x - center[0])**2 + (y - center[1])**2 + (z - center[2])**2)
        
        radii = radii.astype(int)
    else:
        ps = np.abs( np.fft.rfftn(vol) )
        if not return_indices:
            x, y, z = np.indices(ps.shape)
            radii = np.sqrt(x**2 + y**2 + z**2)
            radii = radii.astype(int)
        else:
            [x, y, z] = np.mgrid[-dim[0]//2+m[0]:(dim[0]-1)//2+1, -dim[1]//2+m[1]:(dim[1]-1)//2+1, 0:dim[2]//2+1]
            x = np.fft.ifftshift(x)
            y = np.fft.ifftshift(y)
            radii = np.sqrt(x**2 + y**2 + z**2)
            radii = radii.astype(int)
    radial_profile = np.bincount(radii.ravel(), ps.ravel()) / np.bincount(radii.ravel())
    # exclude corner frequencies
    radial_profile = radial_profile[0:int(round((ps.shape[0]/2)))]
    if not return_indices:
        return radial_profile
    else:
        return radial_profile, radii

def compute_scale_factors(em_profile, ref_profile, apix, scale_factor_arguments, 
                          use_theoretical_profile=True, check_scaling=False):
    
    from locscale.include.emmer.ndimage.profile_tools import scale_profiles, merge_two_profiles
    #print("checkScaling", check_scaling)
    #print("useTheoretical", use_theoretical_profile)
    if use_theoretical_profile:
        theoretical_profile_tuple = get_theoretical_profile(length=len(ref_profile),apix=apix)
        freq = theoretical_profile_tuple[0]
        reference_profile_tuple = (freq, ref_profile)
        
        scaled_theoretical_tuple = scale_profiles(reference_profile_tuple, theoretical_profile_tuple,
                                                  wilson_cutoff=scale_factor_arguments['high_freq'], fsc_cutoff=scale_factor_arguments['fsc_cutoff'])
        
        scaled_theoretical_amplitude = scaled_theoretical_tuple[1]
        smooth = scale_factor_arguments['smooth']
        scaled_reference_profile = merge_two_profiles(ref_profile,scaled_theoretical_amplitude,freq,smooth=smooth,d_cutoff=scale_factor_arguments['wilson'])
        
        reference_profile_for_scaling = scaled_reference_profile
        if check_scaling:
            #print(check_scaling)
            temporary_dictionary = {}
            temporary_dictionary['em_profile'] = em_profile
            temporary_dictionary['input_ref_profile'] = ref_profile
            temporary_dictionary['freq'] = freq
            temporary_dictionary['theoretical_amplitude'] = theoretical_profile_tuple[1]
            temporary_dictionary['scaled_theoretical_amplitude'] = scaled_theoretical_amplitude
            temporary_dictionary['scaled_reference_profile'] = scaled_reference_profile
            temporary_dictionary['scaling_condition'] = scale_factor_arguments
    else:
        reference_profile_for_scaling = ref_profile
        
    np.seterr(divide='ignore', invalid='ignore');
    scale_factor = np.divide(np.abs(reference_profile_for_scaling), np.abs(em_profile))
    scale_factor[ ~ np.isfinite( scale_factor )] = 0; #handle division by zero    
    
    if check_scaling and use_theoretical_profile:
        #print("checkScalingReport", check_scaling)
        temporary_dictionary['scale_factor'] = scale_factor
        return scale_factor, temporary_dictionary
    else:
        return scale_factor

def set_radial_profile(vol, scale_factor, radii):
    ps = np.fft.rfftn(vol)
    for j,r in enumerate(np.unique(radii)[0:vol.shape[0]//2]):
            idx = radii == r
            ps[idx] *= scale_factor[j]

    return np.fft.irfftn(ps, s=vol.shape)

def get_central_scaled_pixel_vals_after_scaling(emmap, modmap, masked_xyz_locs, wn, apix, use_theoretical_profile,scale_factor_arguments, verbose=False,f_cutoff=None, process_name='LocScale', audit=True):
    from tqdm import tqdm
    from locscale.include.emmer.ndimage.map_tools import compute_real_space_correlation
    from locscale.utils.general import true_percent_probability
    import pickle
    
    
    sharpened_vals = []
    central_pix = int(round(wn / 2.0))
    total = (masked_xyz_locs - wn / 2).shape[0]
    cnt = 1.0
    print("Inside get_central_scaled_pixel_vals", use_theoretical_profile)
    
    mpi=False
    if process_name != 'LocScale':
        mpi=True
        from mpi4py import MPI
        
        comm = MPI.COMM_WORLD
        rank=comm.Get_rank()
        size=comm.Get_size()
        
        pbar = {}
        for n in range(size):
            
            pbar[n] = tqdm(total=len(masked_xyz_locs),desc=process_name)
    else:
        progress_bar=tqdm(total=len(masked_xyz_locs), desc=process_name)
    
    
    if audit:
        profiles_audit = {}
    for k, j, i in masked_xyz_locs - wn / 2:
        
        k,j,i,wn = int(round(k)),int(round(j)),int(round(i)),int(round(wn))
        emmap_wn = emmap[k: k+wn, j: j+wn, i: i+ wn]
        modmap_wn = modmap[k: k+wn, j: j+wn, i: i+ wn]

        em_profile = compute_radial_profile(emmap_wn)
        mod_profile, radii = compute_radial_profile(modmap_wn, return_indices=True)
        
        
        check_scaling=true_percent_probability(1) # Checks scaling operation for 1% of all voxels. 
        
        if check_scaling and use_theoretical_profile:
            scale_factors,report = compute_scale_factors(em_profile, mod_profile,apix=apix,scale_factor_arguments=scale_factor_arguments, use_theoretical_profile=use_theoretical_profile,
check_scaling=check_scaling)
            profiles_audit[(k,j,i)] = report
        else:
            scale_factors = compute_scale_factors(em_profile, mod_profile,apix=apix, scale_factor_arguments=scale_factor_arguments, use_theoretical_profile=use_theoretical_profile,
check_scaling=check_scaling)
        
        map_b_sharpened = set_radial_profile(emmap_wn, scale_factors, radii)
        
        #if verbose:
        #    if cnt%1000 == 0:
        #        print ('  {0} {1:.3} percent complete'.format(process_name, (cnt/total)*100))
        

        sharpened_vals.append(map_b_sharpened[central_pix, central_pix, central_pix])
        
        if mpi:
            pbar[rank].update(1)
        else:
            progress_bar.update(n=1)
        
    if mpi:
        comm.barrier()
    if audit and use_theoretical_profile:
        with open("profiles_audit.pickle","wb") as audit:
            pickle.dump(profiles_audit, audit)
        

    return np.array(sharpened_vals, dtype=np.float32)

def put_scaled_voxels_back_in_original_volume_including_padding(sharpened_vals, masked_indices, map_shape):
    map_scaled = np.zeros(np.prod(map_shape))
    map_scaled[masked_indices] = sharpened_vals
    map_scaled = map_scaled.reshape(map_shape)

    return map_scaled

def run_window_function_including_scaling(emmap, modmap, mask, wn, apix, use_theoretical_profile, scale_factor_arguments, 
                                          verbose=False):
    """
    >>> emmap, modmap, mask = setup_test_data()
    >>> scaled_vol = run_window_function_including_scaling(emmap,modmap,mask,wn=10,apix=1.0)
    >>> np.copy(EMNumPy.em2numpy(scaled_vol))[scaled_vol.get_xsize() / 2][scaled_vol.get_ysize() / 2]
    array([ 0.        ,  0.        ,  0.        ,  0.        ,  0.        ,
            0.        ,  0.        ,  0.        ,  0.        ,  0.        ,
            0.12524424,  0.15562208,  0.18547297,  0.24380369,  0.31203741,
            0.46546721,  0.47914436,  0.31334871,  0.28510684,  0.21345402,
            0.17892323,  0.        ,  0.        ,  0.        ,  0.        ,
            0.        ,  0.        ,  0.        ,  0.        ,  0.        ], dtype=float32)
    """
    from locscale.utils.prepare_inputs import get_xyz_locs_and_indices_after_edge_cropping_and_masking
    
    masked_xyz_locs, masked_indices, map_shape = get_xyz_locs_and_indices_after_edge_cropping_and_masking(mask, wn)

    sharpened_vals = get_central_scaled_pixel_vals_after_scaling(emmap, modmap, masked_xyz_locs, wn, apix, use_theoretical_profile,scale_factor_arguments=scale_factor_arguments,verbose=verbose)

    map_scaled = put_scaled_voxels_back_in_original_volume_including_padding(sharpened_vals, masked_indices, map_shape)

    return map_scaled

def split_sequence_evenly(seq, size):
    """
    >>> split_sequence_evenly(list(range(9)), 4)
    [[0, 1], [2, 3, 4], [5, 6], [7, 8]]
    >>> split_sequence_evenly(list(range(18)), 4)
    [[0, 1, 2, 3, 4], [5, 6, 7, 8], [9, 10, 11, 12, 13], [14, 15, 16, 17]]
    """
    newseq = []
    splitsize = 1.0 / size * len(seq)
    for i in range(size):
        newseq.append(seq[int(round(i * splitsize)):int(round((i + 1) * splitsize))])
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


def run_window_function_including_scaling_mpi(emmap, modmap, mask, wn, apix,use_theoretical_profile,
                                              scale_factor_arguments, verbose=False):
    """
    >>> emmap_name, modmap_name, mask_name = setup_test_data_to_files()
    >>> import subprocess
    >>> n = subprocess.call(mpi_cmd.split())
    >>> scaled_vol = get_image('scaled.mrc')
    >>> np.copy(EMNumPy.em2numpy(scaled_vol))[scaled_vol.get_xsize() / 2][scaled_vol.get_ysize() / 2]
    array([ 0.        ,  0.        ,  0.        ,  0.        ,  0.        ,
            0.        ,  0.        ,  0.        ,  0.        ,  0.        ,
            0.12524424,  0.15562208,  0.18547297,  0.24380369,  0.31203741,
            0.46546721,  0.47914436,  0.31334871,  0.28510684,  0.21345402,
            0.17892323,  0.        ,  0.        ,  0.        ,  0.        ,
            0.        ,  0.        ,  0.        ,  0.        ,  0.        ], dtype=float32)
    >>> n = [os.remove(each_file) for each_file in [emmap_name, modmap_name, mask_name, 'scaled.mrc']]
    """
    from mpi4py import MPI
    from locscale.utils.prepare_inputs import get_xyz_locs_and_indices_after_edge_cropping_and_masking
    
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    size = comm.Get_size()

    if rank == 0:
        masked_xyz_locs, masked_indices, map_shape = \
        get_xyz_locs_and_indices_after_edge_cropping_and_masking(mask, wn)

        zs, ys, xs = masked_xyz_locs.T
        zs = split_sequence_evenly(zs, size)
        ys = split_sequence_evenly(ys, size)
        xs = split_sequence_evenly(xs, size)
    else:
        zs = None
        ys = None
        xs = None

    zs = comm.scatter(zs, root=0)
    ys = comm.scatter(ys, root=0)
    xs = comm.scatter(xs, root=0)

    masked_xyz_locs = np.column_stack((zs, ys, xs))

    process_name = 'LocScale process {0} of {1}'.format(rank + 1, size)

    sharpened_vals = get_central_scaled_pixel_vals_after_scaling(emmap, modmap, masked_xyz_locs, wn, apix,use_theoretical_profile=use_theoretical_profile, scale_factor_arguments=scale_factor_arguments,verbose=verbose,process_name=process_name)
    
    sharpened_vals = comm.gather(sharpened_vals, root=0)

    if rank == 0:
        sharpened_vals = merge_sequence_of_sequences(sharpened_vals)

        map_scaled = put_scaled_voxels_back_in_original_volume_including_padding(np.array(sharpened_vals),
        masked_indices, map_shape)
    else:
        map_scaled = None

    comm.barrier()

    return map_scaled, rank

def write_out_final_volume_window_back_if_required(args, wn, window_bleed_and_pad, LocScaleVol, apix):
    from locscale.utils.prepare_inputs import pad_or_crop_volume
    from locscale.include.emmer.ndimage.map_utils import save_as_mrc
        
    if window_bleed_and_pad:
        map_shape = [(LocScaleVol.shape[0] - wn), (LocScaleVol.shape[1] - wn), (LocScaleVol.shape[2] - wn)]
        LocScaleVol = pad_or_crop_volume(LocScaleVol, (map_shape))

        
    save_as_mrc(map_data=LocScaleVol, output_filename=args.outfile, apix=apix, origin=0, verbose=True)
    
    if args.symmetry != "C1":
        print("Imposing a symmetry condition of {}".format(args.symmetry))
        import locscale.include.emda.emda.emda_methods as em
        sym = em.symmetry_average([args.outfile],[args.ref_resolution],pglist=[args.symmetry])
        locscale_symmetry = args.outfile[:-4]+"_{}_symmetry".format(args.symmetry)
        save_as_mrc(map_data=sym[0], output_filename=locscale_symmetry, apix=apix, origin=0, verbose=False)
        print("Find the location of LocScale with symmetry imposed:  {}".format(locscale_symmetry))
        

    return LocScaleVol

