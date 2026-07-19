# Advanced options 

## LocScale arguments
```
optional arguments:
  -h, --help                                       | show this help message and exit
  -em EMMAP_PATH, --emmap_path EMMAP_PATH          | Path to unsharpened EM map
  -hm HALFMAP_PATHS, --halfmap_paths HALFMAP_PATHS | Paths to first and second halfmaps
  -mm MODEL_MAP, --model_map MODEL_MAP             | Path to model map file
  -mc MODEL_COORDINATES, --model_coordinates       | Path to PDB file
  -ma MASK, --mask MASK                            | Input filename mask
  -o OUTFILE, --outfile OUTFILE                    | Output filename
  -v, --verbose                                    | Verbose output
```

 ```
  --output_report                                  | Print a PDF copy of the report
  --report_filename REPORT_FILENAME                | Filename for storing PDF output and statistics
  -op OUTPUT_PROCESSING_FILES, --output_processing_files OUTPUT_PROCESSING_FILES
                        Path to store processing files
```
```
  -wn WINDOW_SIZE, --window_size WINDOW_SIZE
                        window size in pixels
  -mpi, --mpi           MPI version
  -np NUMBER_PROCESSES, --number_processes NUMBER_PROCESSES
                        Number of processes to use
  -ref_it, --refmac_iterations   | For atomic model refinement: number of refmac
                        iterations
  --ref_resolution            | Resolution target for Refmac refinement
  --apix                      | pixel size in Angstrom
  --add_blur                  | Globally sharpen the target map for REFMAC refinement
  --refmac5_path              | Path to refmac5 executable
  --model_resolution          | Resolution limit for Model Map generation
  --symmetry                  | Impose symmetry condition for output
  --fdr_window_size           | Window size in pixels for FDR thresholding
  --fdr_filter                | Pre-filter for FDR thresholding
  --pseudomodel_method        | For pseudomodel: method
  --total_iterations          | For pseudomodel: total iterations
  --distance                  | For pseudomodel: typical distance between atoms
  --molecular_weight          | Input molecular weight (in kDa)
  --build_ca_only             | For gradient pseudomodel building: use only Ca atoms
                                with interatomic distance 3.8
  --smooth_factor             | Smooth factor for merging profiles
  --boost_secondary_structure | Amplify signal corresponding to secondary structures
  --ignore_profiles           | Ignore average secondary structure profile during
                               local scaling
  --no_reference              | Run locscale without using any reference information
  --ignore_profiles           | Ignore average secondary structure profile during
                                local scaling
  --skip_refine               | Ignore REFMAC refinement
```
