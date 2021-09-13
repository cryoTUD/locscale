# -*- coding: utf-8 -*-
"""
Created on Thu Jul  1 11:35:08 2021
"""

import argparse
msg = "using gemmi to simulate EM map from pdb model file, filtered to given resolution and with exact size based on target map"
parser = argparse.ArgumentParser(description=msg)
required_args = parser.add_argument_group('input')
optional_args = parser.add_argument_group('options')
required_args.add_argument("-mdl", "--Modelin", required=True, help='input model (pdb)')
required_args.add_argument("-res", "--resolution", required=True, help='resolution cutoff')
optional_args.add_argument("-map", "--Mapin", required=False, help='target map to determine voxelsize')
optional_args.add_argument("-T", "--transpose", required=False, help='transpose map (XYZ to ZXY convention)')
optional_args.add_argument("-out","--output", required=False, help='filename for generated map')
args = parser.parse_args()

import gemmi
import mrcfile
import numpy as np

def calculate_fourier_frequencies(im, apix):
    """Return the image frequency for every voxel in Fourierspace
    for n-dimensional images
    """
    per_axis_freq = [np.fft.fftfreq(N) for N in im.shape[:-1]]
    per_axis_freq.append(np.fft.rfftfreq(im.shape[-1]))
    dims = np.meshgrid(*per_axis_freq, indexing='ij', sparse=True)
    fourier_frequencies = np.sqrt(np.sum([dim**2 for dim in dims]))
    fourier_frequencies_angstrom = fourier_frequencies / apix
    return fourier_frequencies_angstrom

def tanh_filter(im_freq, cutoff):
    """Returns filter coefficients for a hyperbolic tangent filter. 
    """
    cutoff_freq = 1/cutoff
    filter_fall_off = 0.1;
    filter_coefficients = 1.0 - (1.0 - 0.5*(np.tanh((np.pi*(im_freq+cutoff_freq)/(2*filter_fall_off*cutoff_freq))) - np.tanh((np.pi*(im_freq-cutoff_freq)/(2*filter_fall_off*cutoff_freq)))));
    return filter_coefficients;

def low_pass_filter(im, cutoff, apix):
    """
    Returns a low-pass filter image from a tanh filter.
    """
    im_freq     = calculate_fourier_frequencies(im, apix=apix)
    im_filter   = tanh_filter(im_freq, cutoff);
    im_fft      = np.fft.rfftn(im)
    im_fft_filtered = im_fft * im_filter
    im_filtered = np.fft.irfftn(im_fft_filtered)
    return im_filtered

def save_as_mrc(map_data,voxelsize,cell,outfilename,origin_atom=None):
    with mrcfile.new(outfilename,overwrite=True) as mrc:
        mrc.set_data(np.float32(map_data))
        mrc.voxel_size = voxelsize
        mrc.header.cella.x = cell.x
        mrc.header.cella.y = cell.y
        mrc.header.cella.z = cell.z
        
        if origin_atom is not None:    
            origin = np.rec.array((0,0,0),dtype=[('x','<f4'),('y','<f4'),('z','<f4')])
            origin.x = origin_atom.pdb_position.x
            origin.y = origin_atom.pdb_position.y
            origin.z = origin_atom.pdb_position.z
            print(origin)
            mrc.header.origin = origin
    mrc.close()

pdb_struct = gemmi.read_pdb(args.Modelin)
if args.Mapin:
    target_map_mrc = mrcfile.open(args.Modelin)
    vsize = target_map_mrc.voxel_size.x
    cell = gemmi.UnitCell(a=target_map_mrc.header.cella.x,
                          b=target_map_mrc.header.cella.y,
                          c=target_map_mrc.header.cella.z,
                          alpha=target_map_mrc.header.cellb.alpha,
                          beta=target_map_mrc.header.cellb.beta,
                          gamma=target_map_mrc.header.cellb.gamma)
    pdb_struct.cell = cell
else:
    vsize = 1
    
# prepare dencalc object
pdb_struct.remove_waters()
pdb_struct.remove_empty_chains()    
dencalc = gemmi.DensityCalculatorE()
dencalc.d_min = vsize * 2
dencalc.rate = 1
dencalc.blur= 0
dencalc.set_grid_cell_and_spacegroup(pdb_struct)  
model = pdb_struct[0]
dencalc.put_model_density_on_grid(model)

emmap = np.array(dencalc.grid)
  
if args.resolution is not None:
   emmap = low_pass_filter(emmap,args.resolution,vsize)
 
voxelsize_array = np.rec.array((vsize,vsize,vsize),dtype=[('x','<f4'),('y','<f4'),('z','<f4')])
if args.Mapin and args.output:
    cell = target_map_mrc.header.cella
    save_as_mrc(emmap,voxelsize_array,cell,args.output)
elif args.output:
    cell = np.rec.array((1,1,1),dtype=[('x','<f4'),('y','<f4'),('z','<f4')])
    save_as_mrc(emmap,voxelsize_array,cell,args.output)
    
    