"""
Script to generate reference model map for LocScale scaling

Uses cctbx libraries - please cite:
Grosse-Kunstleve RW et al. J. Appl. Cryst. 35:126-136 (2002)

Arjen Jakobi, EMBL (2016) 
"""

import argparse
import os
import sys
import warnings

from cctbx import crystal, maptbx, xray
from iotbx import ccp4_map, file_reader
import iotbx.pdb
from libtbx import group_args
from mmtbx import real_space_correlation
import mmtbx.maps.utils
from scitbx.array_family import flex

progname = os.path.basename(sys.argv[0])
revision = filter(str.isdigit, "$Revision: 1 $")  # to be updated by gitlab after every commit
datmod = "$Date: 2017-03-06 22:14:31 +0200 (Mo, 06 Mar 2017) $"  # to be updated by gitlab fter every commit
author = 'authors: Arjen J. Jakobi, EMBL' + '; ' + datmod [8:18]
version = progname + '  0.1' + '  (r' + revision + ';' + datmod [6:18] + ')'

simple_cmd = 'phenix.python prepare_locscale_input.py -mc model_coordinates.pdb -em em_map.mrc'

warnings.simplefilter('ignore', DeprecationWarning)

cmdl_parser = argparse.ArgumentParser(
description='*** Computes reference map from PDB model and generates files for LocScale ***\n\n' + \
'Example usage: \"{0}\". '.format(simple_cmd) + \
'Generates copies of input files and simulated map ending with 4locscale.xxx, e.g. em_map_4locscale.mrc ' + \
', model_coordinates_4locscale.pdb and model_coordinates_4locscale.mrc. \n' + \
'{0} on {1}'.format(author, datmod))

cmdl_parser.add_argument('-mc', '--model_coordinates', required=True, type=argparse.FileType('r'), help='Input filename PDB model')
cmdl_parser.add_argument('-em', '--em_map', required=True, type=argparse.FileType('r'), help='Input filename EM map')
cmdl_parser.add_argument('-ma', '--mask', type=argparse.FileType('r'), help='Input filename mask')
cmdl_parser.add_argument('-p', '--apix', type=float, help='pixel size in Angstrom')
cmdl_parser.add_argument('-dmin', '--resolution', type=float, help='map resolution in Angstrom')
cmdl_parser.add_argument('-b', '--b_factor', type=float, default=None, help='set bfactor in [A^2]')
cmdl_parser.add_argument('-a', '--b_add', type=float, default=None, help='add bfactor in [A^2]')
cmdl_parser.add_argument('-r', '--rms', type=float, default=None, help='perturb model by rms')
cmdl_parser.add_argument('-t', '--table', type=str, default="electron", help='Scattering table [electron, itcc]')
cmdl_parser.add_argument('-o', '--outfile', type=str, default="rscc.dat", help='Output filename for RSCC data')

def generate_output_file_names(map, model, mask):
    map_out = os.path.splitext(map.name)[0] + "_4locscale.mrc"
    model_out = os.path.splitext(model.name)[0] + "_4locscale.pdb"
    perturbed_model_out = os.path.splitext(model.name)[0] + "_pertubed.pdb"
    model_map_out = os.path.splitext(model.name)[0] + "_4locscale.mrc"
    mask_out = os.path.splitext(mask.name)[0] + "_4locscale.mrc"
    return map_out, model_out, model_map_out, mask_out, perturbed_model_out

def get_dmin(dmin, target_map):
    if dmin is not None:
        d_min = dmin
        print "Model map will be computed to " + str(d_min) + " Angstrom\n"
    else:
        pixel_size = estimate_pixel_size_from_unit_cell(target_map)
        d_min = round(2 * pixel_size + 0.002, 4)
        print "Model map will be computed to " + str(d_min) + " Angstrom\n"
    return d_min

def check_for_zero_B_factor(xrs):
    xrs = xrs.expand_to_p1(sites_mod_positive=True)
    bs = xrs.extract_u_iso_or_u_equiv()
    sel_zero = bs < 1.e-3
    n_zeros = sel_zero.count(True)
    if (n_zeros > 0):
        print "Input model contains %d atoms with B=0\n" % n_zeros

def print_map_statistics(input_model, target_map):
    print "Map dimensions:", target_map.data.all()
    print "Map origin   :", target_map.data.origin()
    print "Map unit cell:", target_map.unit_cell_parameters
    print ""

def get_symmetry_from_target_map(target_map):
    symm = crystal.symmetry(
    space_group_symbol="P1",
    unit_cell=target_map.unit_cell_parameters)
    return symm    

def estimate_pixel_size_from_unit_cell(target_map):
    unit_cell = target_map.unit_cell_parameters
    map_grid = target_map.data.all()
    apix_x = unit_cell[0] / map_grid[0]
    apix_y = unit_cell[1] / map_grid[1]
    apix_z = unit_cell[2] / map_grid[2]
    try: 
        assert (round(apix_x, 4) == round(apix_y, 4) == round(apix_z, 4))
        pixel_size = round(apix_x, 4)
    except AssertionError:
        print "Inconsistent pixel size: %g" % apix_x  
    return pixel_size

def determine_shift_from_map_header(target_map):
    origin = target_map.data.origin()
    translation_vector = [0 - target_map.data.origin()[0], 0 - target_map.data.origin()[1], 0 - target_map.data.origin()[2]]
    return translation_vector 

def shift_map_to_zero_origin(target_map, cg, map_out, return_map=False):
    em_data = target_map.data.as_double()
    em_data = em_data.shift_origin()
    ccp4_map.write_ccp4_map(
        file_name=map_out,
        unit_cell=cg.unit_cell(),
        space_group=cg.space_group(),
        map_data=em_data,
        labels=flex.std_string([""]))
    if return_map is True:
        shifted_map = file_reader.any_file(map_out).file_object
        return em_data, shifted_map

def apply_shift_transformation_to_model(input_model, target_map, symm, pixel_size, model_out):
    sg = symm.space_group()
    uc = symm.unit_cell()
    pdb_hierarchy = input_model.construct_hierarchy().deep_copy()
    atoms = pdb_hierarchy.atoms()
    sites_frac = uc.fractionalize(sites_cart=atoms.extract_xyz())
    if pixel_size is None:
        pixel_size = estimate_pixel_size_from_unit_cell(target_map)
    translation_vector = determine_shift_from_map_header(target_map)
    translation_vector[0], translation_vector[1], translation_vector[2] = (translation_vector[0] * pixel_size) / uc.parameters()[0], (translation_vector[1] * pixel_size) / uc.parameters()[1], (translation_vector[2] * pixel_size) / uc.parameters()[2]
    new_sites = sites_frac + translation_vector
    translation_vector[0], translation_vector[1], translation_vector[2] = translation_vector[0] * uc.parameters()[0], translation_vector[1] * uc.parameters()[1], translation_vector[2] * uc.parameters()[2]
    atoms.set_xyz(uc.orthogonalize(sites_frac=new_sites))
    f = open(model_out, "w")
    f.write(pdb_hierarchy.as_pdb_string(crystal_symmetry=symm))
    f.close()

def set_isotropic_b_factor(xrs,b_factor):
    xrs.convert_to_isotropic()
    xrs = xrs.set_b_iso(value=b_factor)
    return xrs 

def add_isotropic_b_factor(xrs, b_add):
    xrs.shift_us(b_shift=b_add)
    return xrs

def convert_to_isotropic_b(xrs):
    xrs.convert_to_isotropic()
    return xrs

def apply_random_shift_to_coordinates(xrs, rms, shifted_model, symm, model_out):
    xrs = xrs.random_shift_sites(max_shift_cart=rms)
    pdb_hierarchy = shifted_model.construct_hierarchy()
    pdb_hierarchy.adopt_xray_structure(xrs)
    f = open(model_out, "w")
    f.write(pdb_hierarchy.as_pdb_string(crystal_symmetry=symm))
    f.close()
    return xrs

def compute_model_map(xrs, target_map, symm, d_min, table, model_map_out):
    xrs.scattering_type_registry(
    d_min=d_min,
    table=table)
    fc = xrs.structure_factors(d_min=d_min).f_calc()
    cg = maptbx.crystal_gridding(
    unit_cell=symm.unit_cell(),
    space_group_info=symm.space_group_info(),
    pre_determined_n_real=target_map.data.all())
    fc_map = fc.fft_map(
    crystal_gridding=cg).apply_sigma_scaling().real_map_unpadded()
    try:
        assert (fc_map.all() == fc_map.focus() == target_map.data.all())
    except AssertionError:
        print "Different dimension of experimental and simulated model map."
    ccp4_map.write_ccp4_map(
        file_name=model_map_out,
        unit_cell=cg.unit_cell(),
        space_group=cg.space_group(),
        map_data=fc_map,
        labels=flex.std_string([""]))
    return cg, fc_map 

def compute_real_space_correlation_simple(fc_map, em_data):
    cc_overall_cell = flex.linear_correlation(x=em_data.as_1d(),
                      y=fc_map.as_1d()).coefficient()
    print "\nOverall real-space correlation (unit cell)   : %g\n" % cc_overall_cell

def prepare_reference_and_experimental_map_for_locscale (args, out=sys.stdout):
    """
    """  
    map_out, model_out, model_map_out, mask_out, perturbed_model_out = generate_output_file_names(args.em_map, args.model_coordinates, args.mask)
    target_map = file_reader.any_file(args.em_map.name).file_object
    input_model = file_reader.any_file(args.model_coordinates.name).file_object
    mask = file_reader.any_file(args.mask.name).file_object
    d_min = get_dmin(args.resolution, target_map)
    sc_table = args.table
  
    print_map_statistics(input_model, target_map)
    symm = get_symmetry_from_target_map(target_map)
    apply_shift_transformation_to_model(input_model, target_map, symm, args.apix, model_out)
    shifted_model = iotbx.pdb.hierarchy.input(file_name=model_out)
    xrs = shifted_model.xray_structure_simple(crystal_symmetry=symm)
    if args.b_factor is not None:
       xrs = set_isotropic_b_factor(xrs,args.b_factor)
    elif args.b_add is not None:
       xrs = add_isotropic_b_factor(xrs,args.b_add)
    else:
       xrs = convert_to_isotropic_b(xrs) 
    if args.rms is not None:
       print "Perturbing atoms by rms = ", args.rms, "Angstrom\n"
       xrs = apply_random_shift_to_coordinates(xrs, args.rms, shifted_model, symm, perturbed_model_out)
           
    check_for_zero_B_factor(xrs) 
    cg, fc_map = compute_model_map(xrs, target_map, symm, d_min, sc_table, model_map_out)
    
    try: 
        assert (fc_map.all() == fc_map.focus() == target_map.data.all())
    except AssertionError:
        print "Different dimension of experimental and simulated model map."
  
    shift_map_to_zero_origin(mask, cg, mask_out)
    em_data, shifted_map = shift_map_to_zero_origin(target_map, cg, map_out, return_map=True)
    
    compute_real_space_correlation_simple(em_data, fc_map)	

if (__name__ == "__main__") :
    args = cmdl_parser.parse_args()
    prepare_reference_and_experimental_map_for_locscale(args)

