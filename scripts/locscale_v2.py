import os
import sys
import argparse
from argparse import RawTextHelpFormatter
from mpi4py import MPI

progname = os.path.basename(sys.argv[0])
datmod = "2021-09-07"  # to be updated by gitlab after every commit
author = '\n\nAuthors: Arjen J. Jakobi, Carsten Sachse, EMBL and Alok Bharadwaj'
version = progname + '  0.2' + '  (;' + datmod+ ')'

simple_cmd = 'python np_locscale.py -em emmap.mrc -mm modmap.mrc -ma mask.mrc -w 25 -o loc_scaled.mrc'

cmdl_parser = argparse.ArgumentParser(
description='*** Computes contrast-enhanced cryo-EM maps by local amplitude scaling using a reference model ***\n' + \
('\nExample usage: \"{0}\". {1} on {2}'.format(simple_cmd, author, datmod)),formatter_class=RawTextHelpFormatter)


mpi_cmd = 'mpirun -np 4 python locscale_mpi.py -em emmap.mrc -mm modmap.mrc -ma mask.mrc -p 1.0 -w 10 -mpi -o scaled.mrc'

cmdl_parser.add_argument('-em', '--em_map', required=True, help='Input filename EM map')
cmdl_parser.add_argument('-mm', '--model_map', help='Input filename PDB map')
cmdl_parser.add_argument('-p', '--apix', type=float, help='pixel size in Angstrom')
cmdl_parser.add_argument('-ma', '--mask', help='Input filename mask')
cmdl_parser.add_argument('-w', '--window_size', type=int, help='window size in pixel', default=None)
cmdl_parser.add_argument('-res', '--fsc_resolution', required=True, type=float, help='deposited resolution of EM map')
cmdl_parser.add_argument('-o', '--outfile', required=True, help='Output filename')
cmdl_parser.add_argument('-mpi', '--mpi', action='store_true', default=False,
                         help='MPI version call by: \"{0}\"'.format(mpi_cmd))
cmdl_parser.add_argument('-fdr_w', '--fdr_window_size', type=int, help='window size in pixel for FDR thresholding', default=None)
cmdl_parser.add_argument('-dst', '--distance', type=float, help='For pseudo-atomic model: typical distance between atoms', default=1.2)
cmdl_parser.add_argument('-fdr_f', '--fdr_filter', type=float, help='FDR filter for low pass filtering', default=None)
cmdl_parser.add_argument('-pm', '--pseudomodel_method', help='For pseudo-atomic model: method', default='gradient')
cmdl_parser.add_argument('-it_ref', '--refmac_iterations', help='For pseudo-atomic model: number of refmac iterations', default=5)
cmdl_parser.add_argument('-it', '--total_iterations', type=int, help='For pseudo-atomic model: total iterations', default=None)
cmdl_parser.add_argument('-b_global', '--global_sharpen', type=int, help='Globally sharpen the map', default=None)
cmdl_parser.add_argument('-v', '--verbose', default=True,
                         help='Verbose output')
cmdl_parser.add_argument('-use_pm', '--use_pseudomaps', default=True,
                         help='Use pseudo-atomic model')
def launch_amplitude_scaling(args):
    
    from utils.preparation import prepare_mask_and_maps_for_scaling
    from utils.compute import run_window_function_including_scaling, run_window_function_including_scaling_mpi
    from utils.compute import write_out_final_volume_window_back_if_required
    
    if args.verbose and not args.mpi:
        print('\n  LocScale Arguments\n')
        for arg in vars(args):
            print('    {} : {}'.format(arg, getattr(args, arg)))
            
    
    if not args.mpi:
        parsed_arguments = prepare_mask_and_maps_for_scaling(args)
        #emmap, modmap, mask, wn, window_bleed_and_pad, apix, use_pseudomaps, wilson_cutoff, fsc_cutoff, verbose = prepare_mask_and_maps_for_scaling(args)
        ## ^ Above order is changed
    elif args.mpi:
        comm = MPI.COMM_WORLD
        rank = comm.Get_rank()
        size = comm.Get_size()
        
        if rank==0:
            parsed_arguments = prepare_mask_and_maps_for_scaling(args)
            #emmap, modmap, mask, wn, window_bleed_and_pad, apix, use_pseudomaps, wilson_cutoff, fsc_cutoff, verbose = prepare_mask_and_maps_for_scaling(args)
            ## ^ Above order is changed
        else:
            parsed_arguments = None
        
        #parsed_arguments = comm.scatter(parsed_arguments, root=0)
        comm.barrier()
        parsed_arguments = comm.gather(parsed_arguments, root=0)
            
        
    
    
    
    use_pseudomaps = bool(args.use_pseudomaps)
    
    if not args.mpi:
        input_to_scaling = parsed_arguments[:-1]
        LocScaleVol = run_window_function_including_scaling(*input_to_scaling)
        #LocScaleVol = run_window_function_including_scaling(emmap, modmap, mask, wn, apix, use_theoretical_profile=use_pseudomaps, 
                                                #            wilson_cutoff=wilson_cutoff, fsc_cutoff=fsc_cutoff, verbose=args.verbose)
        
    elif args.mpi:
        input_to_scaling = parsed_arguments[:-1]
        LocScaleVol, rank = run_window_function_including_scaling_mpi(*input_to_scaling)
        #LocScaleVol, rank = run_window_function_including_scaling_mpi(emmap, modmap, mask, wn, apix, 
        #                                                              use_theoretical_profile=use_pseudomaps,
        #                                                              wilson_cutoff=wilson_cutoff, fsc_cutoff=fsc_cutoff, verbose=args.verbose)
                                                                      
        
    
    window_bleed_and_pad = parsed_arguments[-1]
    wn = parsed_arguments[3]
    apix = parsed_arguments[4]
    write_out_final_volume_window_back_if_required(args, wn, window_bleed_and_pad, LocScaleVol, apix=apix)
    
    print("You can find the scaled map here: {}".format(args.outfile))


def main():
    args = cmdl_parser.parse_args()

    launch_amplitude_scaling(args)

if __name__ == '__main__':
    main()
