'''
This module is an extension to LocScale code (np_locscale_fft.py): 
Includes option -fsc to calculate local FSC between two maps.
'''

import np_locscale_fft
from np_locscale_avg import *
import fft_calc
from fft_calc import pyfftw_flag
from fsc_calc import *
from fsc_attributes_calc import *
import mrcfile
import numpy as np
import math, sys

def compute_radial_profile_fsc(vol, center=[0,0,0],return_indices=False,
                               pyfftwobj=None, bytealigned=None):
    dim = vol.shape
    m = np.mod(vol.shape,2)
    # make compliant with both fftn and rfftn
    if center is None:
        ps = np.fft.fftshift((fft_calc.calculate_fft(vol,keep_shape=True)))
        z, y, x = np.indices(ps.shape)
        center = tuple((a - 1) / 2.0 for a in ps.shape[::-1])
        radii = np.sqrt((x - center[0])**2 + (y - center[1])**2 + (z - center[2])**2)
        radii = radii.astype(np.int)
    else:
        if pyfftw_flag:
            ps = pyfftwobj(bytealigned)
        else:
            ps = fft_calc.calculate_fft(vol)
        if not return_indices:
            x, y, z = np.indices(ps.shape)
            radii = np.sqrt(x**2 + y**2 + z**2)
            radii = radii.astype(np.int)
        else:
            [x, y, z] = np.mgrid[-dim[0]//2+m[0]:(dim[0]-1)//2+1, 
                                 -dim[1]//2+m[1]:(dim[1]-1)//2+1, 
                                 0:dim[2]//2+1]
            x = np.fft.ifftshift(x)
            y = np.fft.ifftshift(y)
            radii = np.sqrt(x**2 + y**2 + z**2)
            radii = radii.astype(np.int)
    if return_indices: 
        return ps, radii
    else: return ps

def get_central_pixel_vals_fsc(emmap, modmap, masked_xyz_locs, wn, apix,
                                                verbose=False, 
                                                process_name='LocScale',
                                                avg=False,
                                                maxRes=1.5,
                                                minRes=20.0,
                                                auc=False,
                                                randomonly=False,
                                                nrand=250,
                                                list_medians=None,
                                                softedge=True):
    '''
    Given two maps, calculate FSCs on sliding windows (size wn)
    Only points within the mask is used for computation
    Edges of windows are smoothed
    A random set of 200 fscs can be calculated (randomonly=True)
    
    '''
    fsc_vals = []
    central_pix = int(round(wn / 2.0))
    total = (masked_xyz_locs - wn / 2).shape[0]
    cnt = 1.0
    if randomonly:
        print "Random sampling local FSCs to set min and max resolution range"
    random_selection = random.sample(range(1,total),nrand)
    dict_random_fscs = {}
    sols=[]
    list_sols = []
    cnt_bad_sols = 0
    list_weights = None
    #For pyfftw fft
    em_fftoutput = mod_fftoutput = None
    em_inputarray = mod_inputarray = None
    #soft edging
    if softedge:
        edge = calculate_edge(wn)
        softedge_filter = make_soft_edged_window((wn,wn,wn),edge=edge)
    emmap_wn = np.zeros((wn,wn,wn),dtype=np.float32)
    modmap_wn = np.zeros((wn,wn,wn),dtype=np.float32)
    for k, j, i in (masked_xyz_locs - wn / 2):
        if randomonly:
            if cnt > 1 and int(cnt) not in random_selection: 
                cnt += 1
                continue
        emmap_wn_slice = emmap[k: k+wn, j: j+wn, i: i+ wn]
        modmap_wn_slice = modmap[k: k+wn, j: j+wn, i: i+ wn]
        
        if not compare_tuple(emmap_wn_slice.shape,softedge_filter.shape):
            softedge_filter = make_soft_edged_window(emmap_wn_slice.shape,edge=edge)
            emmap_wn = np.zeros(emmap_wn_slice.shape,dtype=np.float32)
            modmap_wn = np.zeros(modmap_wn_slice.shape,dtype=np.float32)
            
        emmap_wn[:] = emmap_wn_slice*softedge_filter
        modmap_wn[:] = modmap_wn_slice*softedge_filter
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

        em_ps = compute_radial_profile_fsc(emmap_wn,
                                           pyfftwobj=em_pyfftwobj,
                                           bytealigned=em_inputarray)
        mod_ps, radii = compute_radial_profile_fsc(modmap_wn, 
                                                   return_indices=True,
                                                   pyfftwobj=mod_pyfftwobj,
                                                   bytealigned=mod_inputarray)
        list_freq, list_fsc, list_nsf = calculate_fsc(em_ps, mod_ps, radii, 
                                            emmap_wn.shape)
        if len(list_fsc) == 0: 
            resolution = 0.
        elif all(c == 0.0 for c in list_fsc): 
            resolution = 0.0
        elif auc and not randomonly:
            #calculate weights by comparing current fsc to medians from random set
            list_weights = calculate_weights_from_median_deviation(
                                                            list_fsc, 
                                                            list_medians) 
            #calculate area under the spline fit fsc curve in the resolution range
            resolution = calc_area_under_fsc_curve(list_freq, list_fsc,
                                                   map_apix=apix,
                                                   minRes=minRes,
                                                   maxRes=maxRes,
                                                   list_weights=list_weights)
        else:
            if not randomonly:
                #from the medians calculated based on random set,
                #calculate weight for fsc value in each shell
                #the weights will be used for spline fitting on fsc curve
                list_weights = calculate_weights_from_median_deviation(
                                                            list_fsc, 
                                                            list_medians)
            #calculate resolution at fsc 0.5
            #sols : list of potential frequencies at fsc0.5 
            #TODO: check if sols always has one value? 
            resolution,sols = get_resolution_from_fsc_curve(
                                                list_freq, list_fsc,
                                                   map_apix=apix,
                                                   minRes=minRes,
                                                   maxRes=maxRes,
                                                   list_weights=list_weights)
        if np.isnan(resolution):
            resolution = apix/maxRes
        fsc_vals.append(resolution)
        #plot random fscs
        if int(cnt) in random_selection:
            #save randomly selected freq vs fsc in a dictionary for plotting
            dict_random_fscs[cnt] = [np.array(list_freq), 
                                     np.array(list_fsc)]
            if len(sols) > 2: cnt_bad_sols += 1
            elif len(sols) == 0: cnt_bad_sols += 1
            #save potential freq at fsc0.5 in a list
            list_sols.extend(sols)
        if not randomonly and verbose:
            if cnt%1000 == 0:
                print ('  {0} {1:.3} percent complete'.format(
                    process_name, (cnt/total)*100))
            sys.stdout.flush()
        cnt += 1
    mad = np.median(np.absolute(np.array(fsc_vals)-np.median(fsc_vals)))
    #final steps if it is a random set calculation
    if randomonly:
        minfreq = apix/minRes
        maxfreq = apix/maxRes
        #more than 20% strange FSC curves
        if cnt_bad_sols/nrand > 0.1:
            warnings.warn('''FSC curves cross 0.5 multiple times
                        check the maps for artifacts (e.g. from tight masks)''')
            return fsc_vals, minRes, maxRes
        #select min and max freq based on selected fsc05 frequencies
        minfreq,maxfreq = select_min_and_max_freq(list_sols, list_freq, n=nrand)
        #extend minfreq until background/solvant? artifacts 
        if auc:
            minfreq_set = min(apix/7.5,apix/(maxRes+5.0))
            minfreq,maxfreq = extend_minfreq_by_deviation(dict_random_fscs,
                                                  apix/minfreq_set,
                                                  maxRes,
                                                  apix)
            if minfreq > minfreq_set: minfreq = 0.
        if maxfreq is not None:
            if apix/maxfreq > maxRes:
                maxRes = apix/maxfreq
        
        #get medians and median absolute deviations from fsc values at each shell
        try: list_medians, list_mad = get_medians_fsc(dict_random_fscs)
        except: list_medians = None
        #set minRes
        if not minfreq is None and minfreq != 0.: 
            if auc: minRes = apix/minfreq
            elif apix/minfreq < minRes:
                minRes = apix/minfreq
            
        print "Minimum and Maximum resolution (random sample): ", minRes, maxRes
        plot_fscs(dict_random_fscs,"locfscs.png",xlabel="resolution",ylabel="FSC",
              map_apix=apix,lstyle=False,marker=False,minRes=minRes,maxRes=maxRes)
        #sprint "Mean and standard deviation of local FSCs (random sample)", np.mean(fsc_vals), np.std(fsc_vals)
        return fsc_vals, minRes, maxRes, list_medians
    #print "Median and deviation of local FSCs: ", np.median(fsc_vals), mad
    if not auc:
        print "Mean and standard deviation of local FSCs (0.5)", np.mean(fsc_vals), np.std(fsc_vals)
    else:
        print "Mean and standard deviation of local FSCs (AUC)", np.mean(fsc_vals), np.std(fsc_vals)
    return fsc_vals

def run_window_function_including_fsc(emmap, modmap, mask, wn, apix, 
                                      verbose=False, avg=False,
                                      minRes=20.0,maxRes=1.5,auc=False,
                                      randomonly=False):
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
    masked_xyz_locs, masked_indices, map_shape = np_locscale_fft.get_xyz_locs_and_indices_after_edge_cropping_and_masking(mask, wn)
    
    random_fsc_vals,minRes,maxRes, list_medians = get_central_pixel_vals_fsc(
                                                emmap, modmap, 
                                                masked_xyz_locs, wn, 
                                                apix, verbose=verbose,
                                                avg=avg,
                                                minRes=minRes,
                                                maxRes=maxRes,
                                                auc=auc,
                                                randomonly=True)
    print 'Min and Max resolution bounds: ', minRes, maxRes
        #sys.exit()
    if not randomonly:
        fsc_vals = get_central_pixel_vals_fsc(emmap, modmap, 
                                                       masked_xyz_locs, wn, 
                                                       apix, verbose=verbose,
                                                       avg=avg,
                                                       minRes=minRes,
                                                       maxRes=maxRes,
                                                       auc=auc,
                                                       list_medians=list_medians,
                                                       randomonly=randomonly)
        map_fsc = \
            np_locscale_fft.put_scaled_voxels_back_in_original_volume_including_padding(fsc_vals, masked_indices, map_shape)
        #set background to max
        if not auc: map_fsc[map_fsc==0.] = np.amax(map_fsc)
        else: map_fsc[map_fsc==0.] = np.amin(map_fsc)
        return map_fsc

def run_window_function_including_fsc_mpi(emmap, modmap, mask, wn, apix,
                                              verbose=False,avg=False,
                                              minRes=20.0,
                                              maxRes=1.5,auc=False):
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
    
    random_fsc_vals,minRes,maxRes, list_medians = get_central_pixel_vals_fsc(
                                                emmap, modmap, 
                                                masked_xyz_locs, wn, apix,
                                                verbose=verbose, 
                                                process_name=process_name, 
                                                avg=avg,
                                                minRes=minRes,
                                                maxRes=maxRes,
                                                auc=auc,
                                                randomonly=True)
    print 'Min and Max resolution bounds: ', minRes, maxRes
    fsc_vals = get_central_pixel_vals_fsc(emmap, modmap, 
                                                   masked_xyz_locs, wn, apix,
                                                   verbose=verbose, 
                                                   process_name=process_name, 
                                                   avg=avg,
                                                   minRes=minRes,
                                                   maxRes=maxRes,
                                                   auc=auc,
                                                   list_medians=list_medians,
                                                   randomonly=randomonly)
    fsc_vals = comm.gather(fsc_vals, root=0)

    if rank == 0:
        fsc_vals = np_locscale_fft.merge_sequence_of_sequences(fsc_vals)
        map_fsc = np_locscale_fft.put_scaled_voxels_back_in_original_volume_including_padding(np.array(fsc_vals),
        masked_indices, map_shape)
        #set background to max
        map_fsc[map_fsc==0.] = np.amax(map_fsc)
    else:
        map_fsc = None
 
    comm.barrier()

    return map_fsc, rank

def launch_fsc(args):
    if args.verbose and not args.mpi:
        print('\n  LocScale Arguments\n')
        for arg in vars(args):
            print('    {} : {}'.format(arg, getattr(args, arg)))
    emmap, modmap, mask, wn, window_bleed_and_pad, emmapobj = \
                                prepare_mask_and_maps_for_scaling(args)
    #if args.minres is None: args.minres = 20.0
    if args.maxres is None: args.maxres = np.round(2*args.apix,2)
    if not args.mpi:
        LocScaleVol = run_window_function_including_fsc(
                                        emmap, modmap, mask, wn, 
                                        args.apix, verbose=args.verbose,
                                        avg=args.avg,
                                        maxRes=args.maxres,
                                        minRes=args.minres,
                                        auc=args.auc,
                                        randomonly=args.random)
        if not args.random:
            LocScaleVol = write_out_final_volume_window_back_if_required(
                                        args, wn, window_bleed_and_pad, 
                                        LocScaleVol, outfile=args.outfile,
                                        emmapobj=emmapobj)
    elif args.mpi:
        LocScaleVol, rank = run_window_function_including_fsc_mpi(
                                        emmap, modmap, mask, wn, 
                                        args.apix, verbose=args.verbose,
                                        avg=args.avg,
                                        maxRes=args.maxres,
                                        minRes=args.minres,
                                        auc=args.auc,
                                        randomonly=args.random)
        if rank == 0:
            LocScaleVol = write_out_final_volume_window_back_if_required(
                                        args, wn, window_bleed_and_pad, 
                                        LocScaleVol,outfile=args.outfile,
                                        emmapobj=emmapobj)

def main():
    np_locscale_fft.cmdl_parser.add_argument('-avg','--avg', action='store_true', default=False, 
                         help='Scale to average instead of reference')
    np_locscale_fft.cmdl_parser.add_argument('-fsc','--fsc', action='store_true', default=False, 
                         help='Calculate FSCs in local windows')
    np_locscale_fft.cmdl_parser.add_argument('-auc','--auc', action='store_true', default=False, 
                         help='Calculate AUC of FSCs')
    np_locscale_fft.cmdl_parser.add_argument('-minres','--minres', default=20.0, required=False, 
                         help='Minimum resolution for calculation of FSCs')
    np_locscale_fft.cmdl_parser.add_argument('-maxres','--maxres', default=None, required=False, 
                         help='Maximum resolution for calculation of FSCs')
    np_locscale_fft.cmdl_parser.add_argument('-random','--random', action='store_true', default=False, 
                         help='Only a random set of local FSCs')
    args = np_locscale_fft.cmdl_parser.parse_args()
    if args.mpi: args.random = False
    launch_fsc(args)

if __name__ == '__main__':
    sys.exit(main())
