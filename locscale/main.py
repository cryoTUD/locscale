import os
import sys
from locscale.utils.startup_utils import print_start_banner, print_end_banner, print_arguments
from locscale.utils.parse_utils import main_parser
from datetime import datetime


progname = os.path.basename(sys.argv[0])
author = '\n\nAuthors: Arjen J. Jakobi (TU Delft), Alok Bharadwaj (TU Delft), Reinier de Bruin (TU Delft) \n\n'
version = progname + '  2.0'

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
    

if __name__ == '__main__':
    main()
