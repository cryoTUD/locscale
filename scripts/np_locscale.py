import numpy as np
from scipy.ndimage.morphology import binary_dilation
import mrcfile
#import matplotlib.pyplot as plt
import argparse, math, os, sys
from argparse import RawTextHelpFormatter
from locscale_headers import *
import pickle
from emmer.ndimage.profile_tools import *

progname = os.path.basename(sys.argv[0])
datmod = "2018-03-02"  # to be updated by gitlab after every commit
author = '\n\nAuthors: Arjen J. Jakobi and Carsten Sachse, EMBL'
version = progname + '  0.2' + '  (;' + datmod+ ')'

simple_cmd = 'python locscale_mpi.py -em emmap.mrc -mm modmap.mrc -ma mask.mrc -p 1.0 -w 10 -o scaled.mrc'

cmdl_parser = argparse.ArgumentParser(
description='*** Computes contrast-enhanced cryo-EM maps by local amplitude scaling using a reference model ***\n' + \
('\nExample usage: \"{0}\". {1} on {2}'.format(simple_cmd, author, datmod)),formatter_class=RawTextHelpFormatter)


mpi_cmd = 'mpirun -np 4 python locscale_mpi.py -em emmap.mrc -mm modmap.mrc -ma mask.mrc -p 1.0 -w 10 -mpi -o scaled.mrc'

cmdl_parser.add_argument('-em', '--em_map', required=True, help='Input filename EM map')
cmdl_parser.add_argument('-mm', '--model_map', help='Input filename PDB map')
cmdl_parser.add_argument('-p', '--apix', type=float, required=True, help='pixel size in Angstrom')
cmdl_parser.add_argument('-ma', '--mask', help='Input filename mask')
cmdl_parser.add_argument('-w', '--window_size', type=int, help='window size in pixel')
cmdl_parser.add_argument('-r', '--resolution', type=float, help='deposited resolution of EM map')
cmdl_parser.add_argument('-o', '--outfile', required=True, help='Output filename')
cmdl_parser.add_argument('-mpi', '--mpi', action='store_true', default=False,
                         help='MPI version call by: \"{0}\"'.format(mpi_cmd))
cmdl_parser.add_argument('-fdr_w', '--fdr_window_size', type=int, help='window size in pixel for FDR thresholding', default=16)
cmdl_parser.add_argument('-bl', '--bond_length', type=float, help='For pseudo-atomic model: bond length')
cmdl_parser.add_argument('-pm', '--pseudomodel_method', help='For pseudo-atomic model: method')
cmdl_parser.add_argument('-it', '--total_iterations', type=int, help='For pseudo-atomic model: total iterations')

cmdl_parser.add_argument('-v', '--verbose', default=False,
                         help='Verbose output')

use_theoretical_profile_global = True
if use_theoretical_profile_global:
    print("Using theoretical profile")
check_scaling = True

count_to_check_profiles = 0
print(count_to_check_profiles)
if check_scaling:
    random_profiles = {}
    count_to_check_profiles = 1

print(count_to_check_profiles)

# STILL TO BE DONE: replace test data
#def setup_test_data(voldim=30, size=10):
#    from sparx import model_gauss
#    emmap = model_gauss(size, voldim, voldim, voldim)
#    modmap = EMData()
#    modmap.set_size(voldim, voldim, voldim)
#    modmap.process_inplace("testimage.noise.gauss", {"sigma":1, "seed":99})
#    mask = model_square(size, voldim, voldim, voldim)

#    return emmap, modmap, mask

#def setup_test_data_to_files(emmap_name='emmap.mrc', modmap_name='modmap.mrc', mask_name='mask.mrc'):
#    """
#    >>> emmap_name, modmap_name, mask_name = setup_test_data_to_files()
#    >>> import subprocess
#    >>> n = subprocess.call(simple_cmd.split())
#    >>> scaled_vol = get_image('scaled.mrc')
#    >>> np.copy(EMNumPy.em2numpy(scaled_vol))[scaled_vol.get_xsize() / 2][scaled_vol.get_ysize() / 2]
#    array([ 0.        ,  0.        ,  0.        ,  0.        ,  0.        ,
#            0.        ,  0.        ,  0.        ,  0.        ,  0.        ,
#            0.12524424,  0.15562208,  0.18547297,  0.24380369,  0.31203741,
#            0.46546721,  0.47914436,  0.31334871,  0.28510684,  0.21345402,
#            0.17892323,  0.        ,  0.        ,  0.        ,  0.        ,
#            0.        ,  0.        ,  0.        ,  0.        ,  0.        ], dtype=float32)
#    >>> n = [os.remove(each_file) for each_Average bond length forfile in [emmap_name, modmap_name, mask_name, 'scaled.mrc']]
#    """
#    emmap, modmap, mask = setup_test_data()

#    emmap.write_image(emmap_name)
#    modmap.write_image(modmap_name)
#    mask.write_image(mask_name)

#    return emmap_name, modmap_name, mask_name

def compute_padding_average(vol, mask):
    mask = (mask > 0.5).astype(np.int8)
    #inverted_mask = np.logical_not(mask)
    average_padding_intensity = np.mean(np.ma.masked_array(vol, mask))
    return average_padding_intensity

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

def get_theoretical_profile(length,apix):
    
    current_directory = os.getcwd()
    filename = path_to_locscale + "/scripts/theoretical_profiles.pickle"
    with open(filename,'rb') as f:
        profiles = pickle.load(f)
    
    frequency_limits = (float(1/(apix*length)),float(1/(apix*2)))
    helix_profile = profiles['helix']
    resampled_helix_profile = resample_1d(helix_profile['freq'], helix_profile['profile'],num=length,
                                          xlims=frequency_limits)
    return resampled_helix_profile

def get_modmap_from_pseudomodel(args):
    emmap_mrc = mrcfile.open(args.em_map)
    resolution = float(args.resolution)
    fdr_window_size = int(args.fdr_window_size)
    pam_bond_length = float(args.bond_length)
    pam_method = args.pseudomodel_method
    pam_iteration = int(args.total_iterations)
    verbose = bool(args.verbose)
    if args.mask is None:
        print("You have not entered the mask path. Running FDR automatically, to obtain a mask")
        mask_path = run_FDR(emmap_path=args.em_map, window_size=fdr_window_size,verbose=verbose)
        if mask_path is None:
            print("Problem running FDR. Returning None")
            return None
    else:
        mask_path = args.mask
    
    
    num_atoms,mask_dims = measure_mask_parameters(mask_path,verbose=verbose)
    
    pseudomodel_path = run_pam(emmap_path=args.em_map, mask_path=mask_path, threshold=1, num_atoms=num_atoms, 
                               method=pam_method, bl=pam_bond_length,total_iterations=pam_iteration,verbose=verbose)
    if pseudomodel_path is None:
        print("Problem running pseudo-atomic model generator. Returning None")
        return None
    
    refined_model_path = run_refmac(model_path=pseudomodel_path, model_name=pseudomodel_path[:-4], 
                                    map_path=args.em_map, resolution=resolution, maskdims=mask_dims,verbose=verbose)
    if refined_model_path is None:
        print("Problem running REFMAC. Returning None")
        return None
    
    #emmap_path, mask_path = run_mapmask(args.em_map), run_mapmask(mask_path)
    pseudomodel_modmap,new_emmap_path,new_mask_path = run_refmap(model_path=refined_model_path, 
                                                                 emmap_path=args.em_map, 
                                                                 mask_path=mask_path, verbose=verbose)
    
    if pseudomodel_modmap is None:
        print("Problem simulating map from refined model. Returning None")
        return None
    
    # MAPMASK
    pseudomodel_modmap_xyz = run_mapmask(pseudomodel_modmap)
    new_emmap_path_xyz = run_mapmask(new_emmap_path)
    new_mask_path_xyz = run_mapmask(new_mask_path)
    
    
    if pseudomodel_modmap_xyz is None:
        print("Could not finish MAPMASK")
        return None
    else:
        print("Successfully created model map")
        
        return pseudomodel_modmap_xyz,new_emmap_path_xyz,new_mask_path_xyz

    
    
    
def prepare_mask_and_maps_for_scaling(args):
    
    ## Added new part to include pseudo-atomic model
    mask_path = None
    if args.model_map is not None:
        modmap = mrcfile.open(args.model_map).data
        emmap = mrcfile.open(args.em_map).data
    else:
    
        print("Reference Data not supplied! Using pseudo-model")
        modmap_path,emmap_path,mask_path = get_modmap_from_pseudomodel(args)
        emmap = mrcfile.open(emmap_path).data
        mask = (mrcfile.open(mask_path).data>0.5).astype(np.int8)
        modmap =mrcfile.open(modmap_path).data
    #print(modmap.shape)
    
    if args.mask is None:
        mask = np.zeros(emmap.shape)

        if mask.shape[0] == mask.shape[1] and mask.shape[0] == mask.shape[2] and mask.shape[1] == mask.shape[2]:
            rad = mask.shape[0] // 2
            z,y,x = np.ogrid[-rad: rad+1, -rad: rad+1, -rad: rad+1]
            mask = (x**2+y**2+z**2 <= rad**2).astype(np.int_).astype(np.int8)
            mask = pad_or_crop_volume(mask,emmap.shape)
            mask = (mask > 0.5).astype(np.int8)
        else:
            mask += 1
            mask = mask[0:mask.shape[0]-1, 0:mask.shape[1]-1, 0:mask.shape[2]-1]
            mask = pad_or_crop_volume(emmap, (emmap.shape), pad_value=0)
    elif args.mask is not None:
        mask = (mrcfile.open(args.mask).data > 0.5).astype(np.int8)

    #mask = low_pass_filter(mask, 10, mrcfile.open(args.em_map).voxel_size.x)
    
    if args.window_size is None:
        wn = int(math.ceil(round((7 * 3 * args.apix)) /2.) * 2)
    elif args.window_size is not None:
        wn = int(math.ceil(args.window_size / 2.) * 2)

    window_bleed_and_pad = check_for_window_bleeding(mask, wn)
    if window_bleed_and_pad:
        pad_int_emmap = compute_padding_average(emmap, mask)
        pad_int_modmap = compute_padding_average(modmap, mask)
        map_shape = [(emmap.shape[0] + wn), (emmap.shape[1] + wn), (emmap.shape[2] + wn)]
        emmap = pad_or_crop_volume(emmap, map_shape, pad_int_emmap)
        modmap = pad_or_crop_volume(modmap, map_shape, pad_int_modmap)
        mask = pad_or_crop_volume(mask, map_shape, 0)

    return emmap, modmap, mask, wn, window_bleed_and_pad

def compute_radial_profile(vol, center=[0,0,0], return_indices=False):
    dim = vol.shape
    m = np.mod(vol.shape,2)
    # make compliant with both fftn and rfftn
    if center is None:
        ps = np.abs(np.fft.fftshift((np.fft.fftn(vol))))
        z, y, x = np.indices(ps.shape)
        center = tuple((a - 1) / 2.0 for a in ps.shape[::-1])
        radii = np.sqrt((x - center[0])**2 + (y - center[1])**2 + (z - center[2])**2)
        radii = radii.astype(np.int)
    else:
        ps = np.abs( np.fft.rfftn(vol) )
        if not return_indices:
            x, y, z = np.indices(ps.shape)
            radii = np.sqrt(x**2 + y**2 + z**2)
            radii = radii.astype(np.int)
        else:
            [x, y, z] = np.mgrid[-dim[0]//2+m[0]:(dim[0]-1)//2+1, -dim[1]//2+m[1]:(dim[1]-1)//2+1, 0:dim[2]//2+1]
            x = np.fft.ifftshift(x)
            y = np.fft.ifftshift(y)
            radii = np.sqrt(x**2 + y**2 + z**2)
            radii = radii.astype(np.int)
    radial_profile = np.bincount(radii.ravel(), ps.ravel()) / np.bincount(radii.ravel())
    # exclude corner frequencies
    radial_profile = radial_profile[0:int(round((ps.shape[0]/2)))]
    if not return_indices:
        return radial_profile
    else:
        return radial_profile, radii

def compute_scale_factors(em_profile, ref_profile, f_cutoff=None, apix=None, use_theoretical_profile=False,pos=None):
    if check_scaling:
        temp_dictionary = {}
        temp_dictionary['em_profile'] = em_profile
        temp_dictionary['reference_profile'] = ref_profile
    if use_theoretical_profile:
        theoretical_profile = get_theoretical_profile(length=len(ref_profile),apix=apix)
        scaled_theoretical = scale_profiles((theoretical_profile[0],ref_profile),theoretical_profile,using_reference_profile=False)[1]
        if f_cutoff is None:
                f_cutoff = 0.15
        ref_profile = merge_two_profiles(ref_profile,scaled_theoretical,theoretical_profile[0],smooth=0.3,f_cutoff=f_cutoff)
    if check_scaling:
        temp_dictionary['scaled_profile'] = ref_profile
        temp_dictionary['scaled_theoretical'] = scaled_theoretical
        random_profiles[pos] = temp_dictionary
    
    scale_factor = np.divide(np.abs(ref_profile), np.abs(em_profile))
    scale_factor[ ~ np.isfinite( scale_factor )] = 0; #handle division by zero    
    #scale_factor = np.sqrt(ref_profile**2/em_profile**2)
    return scale_factor

def set_radial_profile(vol, scale_factor, radii):
    ps = np.fft.rfftn(vol)
    for j,r in enumerate(np.unique(radii)[0:vol.shape[0]//2]):
            idx = radii == r
            ps[idx] *= scale_factor[j]

    return np.fft.irfftn(ps, s=vol.shape)

def get_central_scaled_pixel_vals_after_scaling(emmap, modmap, masked_xyz_locs, wn, apix,use_theoretical_profile,
                                                verbose=False,f_cutoff=None, process_name='LocScale'):
    sharpened_vals = []
    central_pix = int(round(wn / 2.0))
    total = (masked_xyz_locs - wn / 2).shape[0]
    cnt = 1.0
    for k, j, i in (masked_xyz_locs - wn / 2):
        k,j,i,wn = int(round(k)),int(round(j)),int(round(i)),int(round(wn))
        emmap_wn = emmap[k: k+wn, j: j+wn, i: i+ wn]
        modmap_wn = modmap[k: k+wn, j: j+wn, i: i+ wn]

        em_profile = compute_radial_profile(emmap_wn)
        mod_profile, radii = compute_radial_profile(modmap_wn, return_indices=True)
        
        scale_factors = compute_scale_factors(em_profile, mod_profile,apix=apix,use_theoretical_profile=use_theoretical_profile_global,f_cutoff=f_cutoff,pos=(k,j,i))
        map_b_sharpened = set_radial_profile(emmap_wn, scale_factors, radii)
        
        if verbose:
            if cnt%1000 == 0:
                print ('  {0} {1:.3} percent complete'.format(process_name, (cnt/total)*100))
            cnt += 1

        sharpened_vals.append(map_b_sharpened[central_pix, central_pix, central_pix])

    return np.array(sharpened_vals, dtype=np.float32)

def put_scaled_voxels_back_in_original_volume_including_padding(sharpened_vals, masked_indices, map_shape):
    map_scaled = np.zeros(np.prod(map_shape))
    map_scaled[masked_indices] = sharpened_vals
    map_scaled = map_scaled.reshape(map_shape)

    return map_scaled

def run_window_function_including_scaling(emmap, modmap, mask, wn, apix, use_theoretical_profile,verbose=False,f_cutoff=None):
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
    masked_xyz_locs, masked_indices, map_shape = get_xyz_locs_and_indices_after_edge_cropping_and_masking(mask, wn)

    sharpened_vals = get_central_scaled_pixel_vals_after_scaling(emmap, modmap, masked_xyz_locs, wn, apix, use_theoretical_profile,f_cutoff=f_cutoff,verbose=verbose)

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


def run_window_function_including_scaling_mpi(emmap, modmap, mask, wn, apix,use_theoretical_profile,f_cutoff=None,
                                              verbose=False):
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

    sharpened_vals = get_central_scaled_pixel_vals_after_scaling(emmap, modmap, masked_xyz_locs, wn, apix,
                                                                 verbose=verbose,use_theoretical_profile=use_theoretical_profile_global,f_cutoff=f_cutoff, process_name=process_name)
    sharpened_vals = comm.gather(sharpened_vals, root=0)

    if rank == 0:
        sharpened_vals = merge_sequence_of_sequences(sharpened_vals)

        map_scaled = put_scaled_voxels_back_in_original_volume_including_padding(np.array(sharpened_vals),
        masked_indices, map_shape)
    else:
        map_scaled = None

    comm.barrier()

    return map_scaled, rank

def write_out_final_volume_window_back_if_required(args, wn, window_bleed_and_pad, LocScaleVol):

    if window_bleed_and_pad:
        map_shape = [(LocScaleVol.shape[0] - wn), (LocScaleVol.shape[1] - wn), (LocScaleVol.shape[2] - wn)]
        LocScaleVol = pad_or_crop_volume(LocScaleVol, (map_shape))

    with mrcfile.new(args.outfile) as LocScaleVol_out:
        LocScaleVol_out.set_data(LocScaleVol.astype(np.float32))
        LocScaleVol_out.voxel_size = np.rec.array(( args.apix,  args.apix,  args.apix), dtype=[('x', '<f4'), ('y', '<f4'), ('z', '<f4')])
        LocScaleVol_out.header.nxstart, LocScaleVol_out.header.nystart, LocScaleVol_out.header.nzstart = [0,0,0]

    return LocScaleVol

def launch_amplitude_scaling(args):
    if args.verbose and not args.mpi:
        print('\n  LocScale Arguments\n')
        for arg in vars(args):
            print('    {} : {}'.format(arg, getattr(args, arg)))
    emmap, modmap, mask, wn, window_bleed_and_pad = prepare_mask_and_maps_for_scaling(args)
    if args.model_map is None:
        use_theoretical_profile = True
    else:
        use_theoretical_profile = False
    wilson_cutoff = find_wilson_cutoff(mask_path=args.mask, return_as_frequency=True)

    if not args.mpi:
        LocScaleVol = run_window_function_including_scaling(emmap, modmap, mask, wn, args.apix, use_theoretical_profile=use_theoretical_profile, f_cutoff=wilson_cutoff, verbose=args.verbose)
        LocScaleVol = write_out_final_volume_window_back_if_required(args, wn, window_bleed_and_pad, LocScaleVol)
    elif args.mpi:
        LocScaleVol, rank = run_window_function_including_scaling_mpi(emmap, modmap, mask, wn, args.apix, verbose=args.verbose, use_theoretical_profile=use_theoretical_profile_global,f_cutoff=wilson_cutoff)
        if rank == 0:
            LocScaleVol = write_out_final_volume_window_back_if_required(args, wn, window_bleed_and_pad, LocScaleVol)
    with open('random_profiles_for_check.pickle','wb') as f:
        pickle.dump(random_profiles,f)

def main():
    args = cmdl_parser.parse_args()

    launch_amplitude_scaling(args)

if __name__ == '__main__':
    main()
