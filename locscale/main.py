#
# Delft University of Technology (TU Delft) hereby disclaims all copyright interest in the program 'LocScale'
# written by the Author(s).
# Copyright (C) 2021 Alok Bharadwaj and Arjen J. Jakobi
# This software may be modified and distributed under the terms of the BSD license. 
# You should have received a copy of the BSD 3-clause license along with this program (see LICENSE file file for details).
# If not see https://opensource.org/license/bsd-3-clause/.
#

import os
import sys
from locscale.utils.startup_utils import print_start_banner, print_end_banner, print_arguments
from locscale.utils.parse_utils import main_parser
from datetime import datetime
import locscale

progname = 'locscale'
author = 'Authors: Arjen J. Jakobi (TU Delft), Alok Bharadwaj (TU Delft), Reinier de Bruin (TU Delft)'
version = locscale.__version__

def launch_emmernet(args):
    from locscale.emmernet.utils import check_emmernet_inputs, check_and_save_output
    from locscale.utils.file_tools import change_directory
    from locscale.emmernet.prepare_inputs import prepare_inputs
    from locscale.emmernet.run_emmernet import run_emmernet
    
    ## Print start
    start_time = datetime.now()
    print_start_banner(start_time, "EMmerNet")
    if args.verbose:
        print_arguments(args)
    
    ## Check input
    check_emmernet_inputs(args)
   
    copied_args = change_directory(args, args.output_processing_files)  ## Copy the contents of files into a new directory
    ## Prepare inputs
    input_dictionary = prepare_inputs(copied_args)
    ## Run EMMERNET
    emmernet_output_dictionary = run_emmernet(input_dictionary)
    
    check_and_save_output(input_dictionary, emmernet_output_dictionary)

    ## Print end
    print_end_banner(datetime.now(), start_time)

def launch_locscale(args):
    from locscale.utils.startup_utils import launch_locscale_no_mpi, launch_locscale_mpi
    if args.mpi:
        launch_locscale_mpi(args)
    else:
        launch_locscale_no_mpi(args)
       
        
def test_everything():
    from locscale.tests.utils import download_and_test_everything
    download_and_test_everything()

def main():
    main_args = main_parser.parse_args()
    launch_command = main_args.command

    if launch_command == 'run_emmernet':
        launch_emmernet(main_args)
    elif launch_command == 'run_locscale':
        launch_locscale(main_args)
    elif launch_command == "test":
        test_everything()
    elif launch_command == "version":
        print("LocScale")
        print("Version: ", version)
        print(author)
        print("Python version: ", sys.version)    

if __name__ == '__main__':
    main()
