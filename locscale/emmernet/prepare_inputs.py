import mrcfile
import os

def prepare_inputs(args):
    import os
    from locscale.emmernet.utils import check_emmernet_dependencies, check_and_download_emmernet_model
    from locscale.utils.file_tools import get_emmap_path_from_args
    from locscale.preprocessing.headers import check_axis_order
    print("."*80)

    check_emmernet_dependencies(verbose=True)
    emmernet_model_folder = check_and_download_emmernet_model(verbose=True)

    emmap_path, _ = get_emmap_path_from_args(args)

    ## Get absolute path of emmap_path
    emmap_path_absolute = os.path.abspath(emmap_path)
    ## Get folder of emmap_path
    emmap_folder = os.path.dirname(emmap_path_absolute)

    xyz_emmap_path = check_axis_order(emmap_path)
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
        "emmap_folder": emmap_folder,
        "emmernet_model_folder": emmernet_model_folder,
        }
    
    if verbose:
        print("Inputs parsed successfully")
    print("."*80)
    return inputs_dictionary







