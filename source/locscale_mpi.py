from EMAN2 import EMData, EMNumPy, Util, XYData, Region
import numpy as np
from sparx import get_image, binarize, model_square
import argparse, math, os, sys

progname = os.path.basename(sys.argv[0])
revision = filter(str.isdigit, "$Revision: 1 $")  # to be updated by gitlab after every commit 
datmod = "$Date: 2017-05-30 11:03:31 +0200 (Tu, 30 May 2017) $"  # to be updated by gitlab fter every commit      
author = 'authors: Arjen J. Jakobi and Carsten Sachse, EMBL' + '; ' + datmod [8:18]        
version = progname + '  0.1' + '  (r' + revision + ';' + datmod [7:18] + ')'
     
simple_cmd = 'python locscale.py -em emmap.mrc -mm modmap.mrc -ma mask.mrc -p 1.0 -w 10 -o scaled.mrc'

cmdl_parser = argparse.ArgumentParser(
description='*** Computes contrast-enhanced cryo-EM maps by local amplitude scaling using a reference model ***' + \
'Example usage: \"{0}\". {1} on {2}'.format(simple_cmd, author, datmod))
     
mpi_cmd = 'mpirun -np 4 python locscale.py -em emmap.mrc -mm modmap.mrc -ma mask.mrc -p 1.0 -w 10 -mpi -o scaled.mrc'

cmdl_parser.add_argument('-em', '--em_map', required=True, help='Input filename EM map')
cmdl_parser.add_argument('-mm', '--model_map', required=True, help='Input filename PDB map')
cmdl_parser.add_argument('-p', '--apix', type=float, required=True, help='pixel size in Angstrom')
cmdl_parser.add_argument('-ma', '--mask', help='Input filename mask')
cmdl_parser.add_argument('-w', '--window_size', type=int, help='window size in pixel')
cmdl_parser.add_argument('-o', '--outfile', required=True, help='Output filename')
cmdl_parser.add_argument('-mpi', '--mpi', action='store_true', default=False,
                         help='MPI version call by: \"{0}\"'.format(mpi_cmd))
     
def setup_test_data(voldim=30, size=10):
    from sparx import model_gauss
    emmap = model_gauss(size, voldim, voldim, voldim)
    modmap = EMData()
    modmap.set_size(voldim, voldim, voldim)
    modmap.process_inplace("testimage.noise.gauss", {"sigma":1, "seed":99})
    mask = model_square(size, voldim, voldim, voldim)
    
    return emmap, modmap, mask

def setup_test_data_to_files(emmap_name='emmap.mrc', modmap_name='modmap.mrc', mask_name='mask.mrc'):
    """
    >>> emmap_name, modmap_name, mask_name = setup_test_data_to_files()
    >>> import subprocess
    >>> n = subprocess.call(simple_cmd.split())
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
    emmap, modmap, mask = setup_test_data()

    emmap.write_image(emmap_name)
    modmap.write_image(modmap_name)
    mask.write_image(mask_name)
    
    return emmap_name, modmap_name, mask_name

def map_kurtosis(map):
# requires map as NumPy array
    num = np.sum((map - mean(map)) ** 4)/ len(map)
    denom = variance(map) ** 2  
    return num / denom

def compute_radial_amplitude_distribution(map, apix):
    data = map.do_fft()
    radial_average = data.calc_radial_dist(map.get_xsize() / 2, 0, 1.0, 0)
#     Util.save_data(0,1.0/(apix*map.get_xsize()),radial_average,'test_data.dat')
    return radial_average
  
def set_zero_origin_and_pixel_size(map, apix):
    map['MRC.nxstart'] = 0
    map['MRC.nystart'] = 0
    map['MRC.nzstart'] = 0
    map.set_attr("apix_x", apix)
    map.set_attr("apix_y", apix)
    map.set_attr("apix_z", apix)
    return map
  
def set_radial_amplitude_distribution(map, amplitude_distribution, apix):
    data = map.do_fft()
    frequency_range = np.arange(0, (1 / (2 * apix)), (1.0 / (apix * map.get_xsize())))
    frequency_range = np.ndarray.tolist(frequency_range[0:len(amplitude_distribution)])
    sf = XYData()
    sf.set_xy_list(frequency_range, amplitude_distribution)
    data.process_inplace("filter.setstrucfac", {"apix":apix, "strucfac":sf})
      
    map = data.do_ift()
    return map
  
def get_xyz_locs_and_indices_after_edge_cropping_and_masking(mask, wn):
    mask = np.copy(EMNumPy.em2numpy(mask))
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

def get_central_scaled_pixel_vals_after_scaling(emmap, modmap, masked_xyz_locs, wn, apix):
    sharpened_vals = np.array([], dtype=np.float32)
    
    central_pix = int(round(wn / 2.0))
    for k, j, i in (masked_xyz_locs - wn / 2):
        reg = Region(i, j, k, wn, wn, wn)

        emmap_wn = emmap.get_clip(reg)
        modmap_wn = modmap.get_clip(reg)
        
        mod_radial_average = compute_radial_amplitude_distribution(modmap_wn, apix)
        map_b_sharpened = set_radial_amplitude_distribution(emmap_wn, mod_radial_average, apix)
          
        sharpened_vals = np.append(sharpened_vals, map_b_sharpened[central_pix, central_pix, central_pix])
          
    return sharpened_vals

def put_scaled_voxels_back_in_original_volume_including_padding(sharpened_vals, masked_indices, map_shape):
    map_scaled = np.zeros(np.prod(map_shape))
    map_scaled[masked_indices] = sharpened_vals
    map_scaled = map_scaled.reshape(map_shape)
    
    map_scaled = EMNumPy.numpy2em(np.copy(map_scaled))
    
    return map_scaled

def run_window_function_including_scaling(emmap, modmap, mask, wn, apix):
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
 
    sharpened_vals = get_central_scaled_pixel_vals_after_scaling(emmap, modmap, masked_xyz_locs, wn, apix)
     
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
    
def compute_padding_average(map, mask):
    volume_stats_outside_mask = Util.infomask(map, mask, False)
    average_padding_intensity = volume_stats_outside_mask[0]
    
    return average_padding_intensity
    
def check_for_window_bleeding(mask,wn):
    masked_xyz_locs, masked_indices, mask_shape = get_xyz_locs_and_indices_after_edge_cropping_and_masking(mask, 0)
    
    zs, ys, xs = masked_xyz_locs.T
    nk, nj, ni = mask_shape

    if xs.min() < wn / 2 or xs.max() > (ni - wn / 2) or \
    ys.min() < wn / 2 or ys.max() > (nj - wn / 2) or \
    zs.min() < wn / 2 or zs.max() > (nk - wn / 2):
        window_bleed = True
    else:
        window_bleed = False
    
    return window_bleed     
    
def prepare_mask_and_maps_for_scaling(args):
    emmap = get_image(args.em_map)
    modmap = get_image(args.model_map)

    if args.mask is None:
        mask = EMData()
        xsize, ysize, zsize = emmap.get_xsize(), emmap.get_ysize(), emmap.get_zsize()
        mask.set_size(xsize, ysize, zsize)
        mask.to_zero()
        if xsize == ysize and xsize == zsize and ysize == zsize:
            sphere_radius = xsize // 2
            mask.process_inplace("testimage.circlesphere", {"radius":sphere_radius})
        else:
            mask += 1
            mask = Util.window(mask, xsize - 1, ysize - 1, zsize -1)
            mask = Util.pad(mask, xsize, ysize, zsize, 0, 0, 0, '0')
    elif args.mask is not None:
        mask = binarize(get_image(args.mask), 0.5)
        
    if args.window_size is None:
        wn = int(math.ceil(round((7 * 3 * args.apix)) /2.) * 2) 
    elif args.window_size is not None:
        wn = int(math.ceil(args.window_size / 2.) * 2)

    window_bleed_and_pad = check_for_window_bleeding(mask, wn)
    if window_bleed_and_pad:
        pad_int_emmap = compute_padding_average(emmap, mask)
        pad_int_modmap = compute_padding_average(modmap, mask)

        map_shape = [(emmap.get_xsize() + wn), (emmap.get_ysize() + wn), (emmap.get_zsize() + wn)] 
        emmap = Util.pad(emmap, map_shape[0], map_shape[1], map_shape[2], 0, 0, 0, 'pad_int_emmap')
        modmap = Util.pad(modmap, map_shape[0], map_shape[1], map_shape[2], 0, 0, 0, 'pad_int_modmap')
        mask = Util.pad(mask, map_shape[0], map_shape[1], map_shape[2], 0, 0, 0, '0')
        
    return emmap, modmap, mask, wn, window_bleed_and_pad

def run_window_function_including_scaling_mpi(emmap, modmap, mask, wn, apix):
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
 
    sharpened_vals = get_central_scaled_pixel_vals_after_scaling(emmap, modmap, masked_xyz_locs, wn, apix)
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
    LocScaleVol = set_zero_origin_and_pixel_size(LocScaleVol, args.apix)
    if window_bleed_and_pad:
        map_shape = [(LocScaleVol.get_xsize() - wn), (LocScaleVol.get_ysize() - wn), (LocScaleVol.get_zsize() - wn)]
        LocScaleVol = Util.window(LocScaleVol, map_shape[0], map_shape[1], map_shape[2])

    LocScaleVol.write_image(args.outfile)

    return LocScaleVol

def launch_amplitude_scaling(args):
    emmap, modmap, mask, wn, window_bleed_and_pad = prepare_mask_and_maps_for_scaling(args) 
    if not args.mpi:
        LocScaleVol = run_window_function_including_scaling(emmap, modmap, mask, wn, args.apix)
        LocScaleVol = write_out_final_volume_window_back_if_required(args, wn, window_bleed_and_pad, LocScaleVol)
    elif args.mpi:
        LocScaleVol, rank = run_window_function_including_scaling_mpi(emmap, modmap, mask, wn, args.apix)
        if rank == 0:
            LocScaleVol = write_out_final_volume_window_back_if_required(args, wn, window_bleed_and_pad, LocScaleVol)

def main():
    args = cmdl_parser.parse_args()
    
    launch_amplitude_scaling(args)
  
if __name__ == '__main__':
    main()
