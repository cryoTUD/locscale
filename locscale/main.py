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
from locscale.utils.startup_utils import launch_feature_enhance, launch_contrast_enhance, run_housekeeping, test_everything, print_version
from locscale.utils.parse_utils import locscale_parser
import locscale

run_housekeeping()
    
def main():
    main_args = locscale_parser.parse_args()

    launch_command = main_args.command
    
    if launch_command == 'feature_enhance':
        launch_feature_enhance(main_args)
    elif launch_command == 'version':
        print_version()
    elif launch_command == 'test':
        test_everything()
    elif launch_command is None:
        launch_contrast_enhance(main_args)
    else:
        raise ValueError("Unknown command: ", launch_command)
        

if __name__ == '__main__':
    main()
