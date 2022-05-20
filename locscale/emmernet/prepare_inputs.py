import mrcfile
import os

def prepare_inputs(args):
    from locscale.emmernet.utils import check_emmernet_dependencies
    from locscale.utils.file_tools import get_emmap_path_from_args
    from locscale.preprocessing.headers import run_mapmask
    print("."*80)
    check_emmernet_dependencies(verbose=True)

    emmap_path, _ = get_emmap_path_from_args(args)

    xyz_emmap_path = run_mapmask(emmap_path)
    emmernet_type = args.trained_model
    stride = args.stride
    verbose = args.verbose
    outputfile = args.outfile
    batch_size = args.batch_size
    gpu_ids = args.gpu_ids



    inputs_dictionary = {
        "emmap_path": emmap_path,
        'xyz_emmap_path': xyz_emmap_path,
        "emmernet_type": emmernet_type,
        "stride": stride,
        "verbose": verbose,
        "outfile": outputfile,
        "batch_size": batch_size,
        "gpu_ids": gpu_ids,
        }
    
    if verbose:
        print("Inputs parsed successfully")
    print("."*80)
    return inputs_dictionary







