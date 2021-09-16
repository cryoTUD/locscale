'''
This module is an extension to LocScale code (np_locscale_fft.py): 
Includes option (-avg) to locally scale maps to the average intensity (of the two maps)
 instead of a reference.
'''

import np_locscale_fft
import fft_calc
from fft_calc import pyfftw_flag
import fsc_calc
import mrcfile
import numpy as np
import math, sys

def prepare_mask_and_maps_for_scaling(args):
    #save mrcfile map obj
    emmapobj = mrcfile.open(args.em_map)
    modmapobj = mrcfile.open(args.model_map)
    emmap = emmapobj.data
    modmap = modmapobj.data

    if args.mask is None:
        mask = np.zeros(emmap.shape)

        if mask.shape[0] == mask.shape[1] and mask.shape[0] == mask.shape[2] and mask.shape[1] == mask.shape[2]:
            rad = mask.shape[0] // 2
            z,y,x = np.ogrid[-rad: rad+1, -rad: rad+1, -rad: rad+1]
            mask = (x**2+y**2+z**2 <= rad**2).astype(np.int_).astype(np.int8)
            mask = np_locscale_fft.pad_or_crop_volume(mask,emmap.shape)
            mask = (mask > 0.5).astype(np.int8)
        else:
            mask += 1
            mask = mask[0:mask.shape[0]-1, 0:mask.shape[1]-1, 0:mask.shape[2]-1]
            mask = np_locscale_fft.pad_or_crop_volume(emmap, (emmap.shape), pad_value=0)
    elif args.mask is not None:
        mask = (mrcfile.open(args.mask).data > 0.5).astype(np.int8)

    if args.window_size is None:
        wn = int(math.ceil(round((7 * 3 * args.apix)) /2.) * 2)
    elif args.window_size is not None:
        wn = int(math.ceil(args.window_size / 2.) * 2)

    window_bleed_and_pad = np_locscale_fft.check_for_window_bleeding(mask, wn)
    if window_bleed_and_pad:
        pad_int_emmap = np_locscale_fft.compute_padding_average(emmap, mask)
        pad_int_modmap = np_locscale_fft.compute_padding_average(modmap, mask)
        map_shape = [(emmap.shape[0] + wn), (emmap.shape[1] + wn), (emmap.shape[2] + wn)]
        emmap = np_locscale_fft.pad_or_crop_volume(emmap, map_shape, pad_int_emmap)
        modmap = np_locscale_fft.pad_or_crop_volume(modmap, map_shape, pad_int_modmap)
        mask = np_locscale_fft.pad_or_crop_volume(mask, map_shape, 0)
    #return map obj as well
    return emmap, modmap, mask, wn, window_bleed_and_pad, emmapobj

def compute_scale_factors(em_profile, ref_profile, avg=False):
    '''
    Scale intensity to a reference or average 
    '''
#     #temporarily raise warnings as errors
#     with warnings.catch_warnings():
#         warnings.filterwarnings('error')
    if not avg: 
        scale_factor = np.sqrt(ref_profile**2/em_profile**2)
#             except RuntimeWarning:
#                 #nan/inf as scale factors?
    else: 
        avg_profile = (ref_profile**2+em_profile**2)/2.0
        scale_factor = np.sqrt(avg_profile/em_profile**2)
    #scale_factor if reference is 0 ?
    for l in xrange(len(scale_factor)):
        if np.isnan(scale_factor[l]) or np.isinf(scale_factor[l]): 
            scale_factor[l] = 1.0
    return scale_factor

def get_central_scaled_pixel_vals_after_scaling(emmap, modmap, masked_xyz_locs, 
                                                wn, apix,verbose=False, 
                                                process_name='LocScale',
                                                avg=False,softedge=True):
    sharpened_vals = []
    #avg: initialize for second map
    sharpened_vals_mod = []
    central_pix = int(round(wn / 2.0))
    total = (masked_xyz_locs - wn / 2).shape[0]
    cnt = 1.0
    #For pyfftw fft
    em_fftoutput = mod_fftoutput = None
    em_inputarray = mod_inputarray = None
    #soft edging
    if softedge:
        edge = fsc_calc.calculate_edge(wn)
        softedge_filter = fsc_calc.make_soft_edged_window((wn,wn,wn),edge=edge)
    for k, j, i in (masked_xyz_locs - wn / 2):
        emmap_wn_view = emmap[k: k+wn, j: j+wn, i: i+ wn]
        modmap_wn_view = modmap[k: k+wn, j: j+wn, i: i+ wn]
        #soft edging
        if softedge:
            emmap_wn = np.zeros(emmap_wn_view.shape,dtype=np.float32)
            modmap_wn = np.zeros(modmap_wn_view.shape,dtype=np.float32)
            if not fsc_calc.compare_tuple(emmap_wn.shape,softedge_filter.shape):
                softedge_filter = fsc_calc.make_soft_edged_window(emmap_wn.shape,edge=edge)
            emmap_wn[:] = emmap_wn_view*softedge_filter
            modmap_wn[:] = modmap_wn_view*softedge_filter
        else:
            emmap_wn = emmap_wn_view
            modmap_wn = modmap_wn_view

        if cnt == 1.0:
            if pyfftw_flag:
                #if center of ps is given [0,0,0]
                #plan fft
                try: 
                    em_pyfftwobj, em_fftoutput, em_inputarray = fft_calc.plan_fft(
                                                            emmap_wn,
                                                            new_inparray=True)
                except: em_pyfftwobj = None
                try: 
                    mod_pyfftwobj, mod_fftoutput, mod_inputarray = fft_calc.plan_fft(
                                                            modmap_wn,
                                                            new_inparray=True)
                    
                except: mod_pyfftwobj = None
            else: 
                em_pyfftwobj = None
                mod_pyfftwobj = None
        #pass bytealigned array for pyfftw
        if pyfftw_flag:
            em_inputarray[:,:,:] = emmap_wn
            mod_inputarray[:,:,:] = modmap_wn
        em_profile = np_locscale_fft.compute_radial_profile(emmap_wn,
                                            pyfftwobj=em_pyfftwobj,
                                            bytealigned=em_inputarray)
        mod_profile, radii = np_locscale_fft.compute_radial_profile(modmap_wn, 
                                                    return_indices=True,
                                                    pyfftwobj=mod_pyfftwobj,
                                                    bytealigned=mod_inputarray)

        scale_factors = compute_scale_factors(em_profile, mod_profile,avg=avg)
        #for pyfftw ifft
        if cnt == 1.0:
            #plan ifft
            if pyfftw_flag and not em_fftoutput is None:
                try: 
                    em_pyifftwobj, em_ifftoutput, em_inp = fft_calc.plan_ifft(
                                                            em_fftoutput,
                                                            output_shape=emmap_wn.shape)
                except: 
                    em_pyifftwobj = None
            else: 
                em_pyifftwobj = None
            ##avg : plan for second map
            if pyfftw_flag and not mod_fftoutput is None:
                try: 
                    mod_pyifftwobj, mod_ifftoutput, mod_inp = fft_calc.plan_ifft(
                                                            mod_fftoutput,
                                                            output_shape=modmap_wn.shape)
                except: 
                    mod_pyifftwobj = None
            else:
                mod_pyifftwobj = None
        map_b_sharpened = np_locscale_fft.set_radial_profile(emmap_wn, scale_factors, radii,
                                             pyfftwobj=em_pyfftwobj,
                                             pyifftwobj=em_pyifftwobj,
                                             bytealigned=em_inputarray)
        #avg: scale second map
        if avg:
            #compute scale factors for the ref map (second map)
            scale_factors_mod = compute_scale_factors(mod_profile, em_profile, avg=avg)
            mod_b_sharpened = np_locscale_fft.set_radial_profile(modmap_wn, scale_factors_mod, radii,
                                                 pyfftwobj=mod_pyfftwobj,
                                                 pyifftwobj=mod_pyifftwobj,
                                                 bytealigned=mod_inputarray)
            sharpened_vals_mod.append(mod_b_sharpened[central_pix, central_pix, central_pix])

        if verbose:
            if cnt%1000 == 0:
                print ('  {0} {1:.3} percent complete'.format(process_name, (cnt/total)*100))
            
            sys.stdout.flush()
        cnt += 1
        
        sharpened_vals.append(map_b_sharpened[central_pix, central_pix, central_pix])
    #avg: return second map as well
    if avg: return np.array(sharpened_vals, dtype=np.float32), \
                    np.array(sharpened_vals_mod, dtype=np.float32)
    return np.array(sharpened_vals, dtype=np.float32)

def run_window_function_including_scaling(emmap, modmap, mask, wn, apix, verbose=False, avg=False):
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
    masked_xyz_locs, masked_indices, map_shape = \
            np_locscale_fft.get_xyz_locs_and_indices_after_edge_cropping_and_masking(mask, wn)

    if not avg: 
        sharpened_vals = get_central_scaled_pixel_vals_after_scaling(emmap, modmap, masked_xyz_locs, wn, apix, verbose=verbose)
        map_scaled = np_locscale_fft.put_scaled_voxels_back_in_original_volume_including_padding(sharpened_vals, masked_indices, map_shape)
        return map_scaled
    else:
        sharpened_vals, sharpened_vals_mod = get_central_scaled_pixel_vals_after_scaling(emmap, modmap, masked_xyz_locs, wn, apix, verbose=verbose,avg=avg)
        map_scaled = np_locscale_fft.put_scaled_voxels_back_in_original_volume_including_padding(sharpened_vals, masked_indices, map_shape)
        modmap_scaled = np_locscale_fft.put_scaled_voxels_back_in_original_volume_including_padding(sharpened_vals_mod, masked_indices, map_shape)
        return map_scaled, modmap_scaled

def run_window_function_including_scaling_mpi(emmap, modmap, mask, wn, apix,
                                              verbose=False,avg=False):
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
        np_locscale_fft.get_xyz_locs_and_indices_after_edge_cropping_and_masking(mask, wn)

        zs, ys, xs = masked_xyz_locs.T
        zs = np_locscale_fft.split_sequence_evenly(zs, size)
        ys = np_locscale_fft.split_sequence_evenly(ys, size)
        xs = np_locscale_fft.split_sequence_evenly(xs, size)
    else:
        zs = None
        ys = None
        xs = None

    zs = comm.scatter(zs, root=0)
    ys = comm.scatter(ys, root=0)
    xs = comm.scatter(xs, root=0)

    masked_xyz_locs = np.column_stack((zs, ys, xs))

    process_name = 'LocScale process {0} of {1}'.format(rank + 1, size)

    if not avg: 
        sharpened_vals = get_central_scaled_pixel_vals_after_scaling(emmap, modmap, masked_xyz_locs, wn, apix,
                                                                     verbose=verbose, process_name=process_name)
        sharpened_vals = comm.gather(sharpened_vals, root=0)
    #avg: get scaled values for both maps
    else: 
        sharpened_vals, sharpened_vals_mod = get_central_scaled_pixel_vals_after_scaling(emmap, modmap, masked_xyz_locs, wn, apix,
                                                                 verbose=verbose, process_name=process_name, avg=avg)
        sharpened_vals = comm.gather(sharpened_vals, root=0)
        sharpened_vals_mod = comm.gather(sharpened_vals_mod, root=0)

    if rank == 0:
        sharpened_vals = np_locscale_fft.merge_sequence_of_sequences(sharpened_vals)
        map_scaled = np_locscale_fft.put_scaled_voxels_back_in_original_volume_including_padding(np.array(sharpened_vals),
        masked_indices, map_shape)
        #avg: merge and box second map scaled
        if avg: 
            sharpened_vals_mod = np_locscale_fft.merge_sequence_of_sequences(sharpened_vals_mod)
            modmap_scaled = np_locscale_fft.put_scaled_voxels_back_in_original_volume_including_padding(
                                            np.array(sharpened_vals_mod),
                                            masked_indices, map_shape)
    else:
        map_scaled = None
        modmap_scaled = None

    comm.barrier()
    #avg: return both maps scaled
    if avg: return map_scaled, modmap_scaled, rank
    
    return map_scaled, rank

def write_out_final_volume_window_back_if_required(args, wn, window_bleed_and_pad, LocScaleVol,avg=False,outfile=None,
                                                   emmapobj=None):

    if window_bleed_and_pad:
        map_shape = [(LocScaleVol.shape[0] - wn), (LocScaleVol.shape[1] - wn), (LocScaleVol.shape[2] - wn)]
        LocScaleVol = np_locscale_fft.pad_or_crop_volume(LocScaleVol, (map_shape))

    LocScaleVol_out = mrcfile.new(outfile,overwrite=True)
    LocScaleVol_out.set_data(LocScaleVol.astype(np.float32))
    LocScaleVol_out.voxel_size = np.rec.array(( args.apix,  args.apix,  args.apix), dtype=[('x', '<f4'), ('y', '<f4'), ('z', '<f4')])
    
    if not emmapobj is None:
        LocScaleVol_out.header.origin.x = emmapobj.header.origin.x
        LocScaleVol_out.header.origin.y = emmapobj.header.origin.y
        LocScaleVol_out.header.origin.z = emmapobj.header.origin.z
        LocScaleVol_out.header.nxstart = emmapobj.header.nxstart
        LocScaleVol_out.header.nystart = emmapobj.header.nystart
        LocScaleVol_out.header.nzstart = emmapobj.header.nzstart
    else: LocScaleVol_out.header.nxstart, LocScaleVol_out.header.nystart, LocScaleVol_out.header.nzstart = [0,0,0]
    
    return LocScaleVol_out


def launch_amplitude_scaling(args):
#     if args.verbose and not args.mpi:
#         print('\n  LocScale Arguments\n')
#         for arg in vars(args):
#             print('    {} : {}'.format(arg, getattr(args, arg)))
    emmap, modmap, mask, wn, window_bleed_and_pad, emmapobj = \
                            prepare_mask_and_maps_for_scaling(args)
    if not args.mpi:
        LocScaleVol = run_window_function_including_scaling(
                                        emmap, modmap, mask, wn, 
                                        args.apix, verbose=args.verbose,
                                        avg=args.avg)
        if args.avg:
            LocScaleVol_map = write_out_final_volume_window_back_if_required(
                                        args, wn, window_bleed_and_pad, 
                                        LocScaleVol[0],avg=args.avg,
                                        outfile='.'.join(args.outfile.split('.')[:-1])+'_map1.mrc'
                                        ,emmapobj=emmapobj)
            LocScaleVol_mod = write_out_final_volume_window_back_if_required(
                                        args, wn, window_bleed_and_pad, 
                                        LocScaleVol[1],avg=args.avg,
                                        outfile='.'.join(args.outfile.split('.')[:-1])+'_map2.mrc'
                                        ,emmapobj=emmapobj)
            return LocScaleVol_map, LocScaleVol_mod
        else:
            LocScaleVol = write_out_final_volume_window_back_if_required(
                                        args, wn, window_bleed_and_pad, 
                                        LocScaleVol, outfile=args.outfile
                                        ,emmapobj=emmapobj)
            return LocScaleVol
    elif args.mpi:
        if args.avg:
            LocScaleVol1, LocScaleVol2,rank = run_window_function_including_scaling_mpi(
                                        emmap, modmap, mask, wn, 
                                        args.apix, verbose=args.verbose,
                                        avg=args.avg)
            if rank == 0:
                LocScaleVol_map = write_out_final_volume_window_back_if_required(
                                        args, wn, window_bleed_and_pad, 
                                        LocScaleVol1,avg=args.avg,
                                        outfile='.'.join(args.outfile.split('.')[:-1])+'_map1.mrc'
                                        ,emmapobj=emmapobj)
                LocScaleVol_mod = write_out_final_volume_window_back_if_required(
                                        args, wn, window_bleed_and_pad, 
                                        LocScaleVol2,avg=args.avg,
                                        outfile='.'.join(args.outfile.split('.')[:-1])+'_map2.mrc'
                                        ,emmapobj=emmapobj)
                return LocScaleVol_map, LocScaleVol_mod
        else:
            LocScaleVol,rank = run_window_function_including_scaling_mpi(
                                        emmap, modmap, mask, wn, 
                                        args.apix, verbose=args.verbose,
                                        avg=args.avg)
            if rank == 0:
                LocScaleVol = write_out_final_volume_window_back_if_required(
                                        args, wn, window_bleed_and_pad, 
                                        LocScaleVol,outfile=args.outfile
                                        ,emmapobj=emmapobj)
                return LocScaleVol

def main():
    np_locscale_fft.cmdl_parser.add_argument('-avg','--avg', action='store_true', default=False, 
                         help='Scale to average instead of reference')
    args = np_locscale_fft.cmdl_parser.parse_args()
    launch_amplitude_scaling(args)

if __name__ == '__main__':
    sys.exit(main())
