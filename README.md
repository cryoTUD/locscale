[![stability-beta](https://img.shields.io/badge/stability-beta-33bbff.svg)](https://github.com/mkenney/software-guides/blob/master/STABILITY-BADGES.md#beta)
[![Python 3.6](https://img.shields.io/badge/python-3.7%20%7C%203.8-brightgreen)](https://www.python.org/downloads/release/python-370/)
[![PyPI](https://img.shields.io/pypi/v/locscale.svg?style=flat)](https://pypi.org/project/locscale/)
[![PyPI - Downloads](https://img.shields.io/pypi/dm/locscale)](https://pypi.org/project/locscale/)
[![License](https://img.shields.io/pypi/l/locscale.svg?color=orange)](https://gitlab.tudelft.nl/aj-lab/locscale/raw/master/LICENSE)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.6652013.svg)](https://doi.org/10.5281/zenodo.6652013)
[![Citation Badge](https://api.juleskreuer.eu/citation-badge.php?doi=10.7554/eLife.27131)](https://juleskreuer.eu/projekte/citation-badge/)

# LocScale - reference-based local sharpening of cryo-EM maps

`LocScale` is an automated program for local sharpening of cryo-EM maps with the aim to improve their interpretability. It utilises general properties inherent to electron scattering from biological macromolecules to restrain the sharpening filter. These can be provided either from an existing atomic model, or inferred directly from the experimental density map.

#### New in LocScale 2.2.3:
- Completely automated process for local map sharpening 

- [Feature_enhance](#4-confidence-aware-density-modification): a confidence-aware density modification tool to enhance features in cryo-EM maps using the `EMmerNet` neural network.

- [Hybrid sharpening](#2-run-locscale-using-a-partial-atomic-model): `LocScale` now supports reference-based sharpening when only partial atomic model information is available.

- [Model-free sharpening](#3-run-locscale-without-atomic-model): `LocScale` now supports reference-based sharpening without the need to supply any atomic model information.

- Full support for point group symmetry (helical symmetry to follow).

<br>
  
`LocScale` is distributed as a portable stand-alone installation that includes all the needed libraries from: https://gitlab.tudelft.nl/aj-lab/locscale/releases.   


Please note that there is a GUI implemented version available as part of the [CCP-EM](http://www.ccpem.ac.uk/) project; LocScale2 is also implemented in [Scipion](http://scipion.i2pc.es/)(thanks to Grigory Sharov, MRC-LMB). Note that currently the CCPEM GUI implementations only supports an older version of Locscale (Locscale 1.0, with only model-based sharpening). 

## Installation 

We recommend to use [Conda](https://docs.conda.io/en/latest/) for a local working environment. See [here](https://docs.conda.io/projects/conda/en/latest/user-guide/install/download.html#anaconda-or-miniconda) for more information on what Conda flavour may be the right choice for you, and [here](https://www.anaconda.com/products/distribution) for Conda installation instructions.

#### Requirements

LocScale should run on any CPU system with Linux, OS X or Windows subsytem for Linux (WSL). To run LocScale efficiently in EMmerNet mode requires the availability of a GPU; it is possible to run it on CPUs but computation will be slow(er). 

#### Installation instructions:

#### Quick installation: 
##### 1. Install REFMAC5 via CCP4/CCPEM
LocScale needs a working instance of [REFMAC5](https://www2.mrc-lmb.cam.ac.uk/groups/murshudov/index.html). If you already have CCP4/CCPEM installed check if the path to run `refmac5` is present in your environment. 

```bash
which refmac5
```

If no valid path is returned, please install [CCP4](https://www.ccp4.ac.uk/download/) to ensure refmac5 is accessible to the program. 

##### 2. Install LocScale using environment files 

There are two yml files in the repo: environment_cpu.yml and environment_gpu.yml. We recommend you to download and install the GPU version.

Once you download the yml file of your choice: 
```bash
conda env create -f /path/to/environment_cpu.yml
conda activate cpu_locscale
```
or 
```bash
conda env create -f /path/to/environment_gpu.yml
conda activate gpu_locscale
```
#### Alternatively
You can also follow these steps to install locscale using pip.

##### 1. Create and activate a new conda environment

```bash
conda create -n locscale python=3.8 
conda activate locscale
```

##### 2. Install fortran compiler
LocScale uses Fortran code to perform symmetry operations and requires a Fortran compiler to be present in your system. You can install `gfortran` from conda-forge.
```bash
conda install -c conda-forge gfortran
```
##### 3. Install REFMAC5 via CCP4/CCPEM

The model-based and hybrid map sharpening modes of LocScale need a working instance of [REFMAC5](https://www2.mrc-lmb.cam.ac.uk/groups/murshudov/index.html). If you already have CCP4/CCPEM installed check if the path to run `refmac5` is present in your environment. For model-free sharpening and confidence-aware density modification REFMAC5 is not required. 

```bash
which refmac5
```

If no valid path is returned, please install [CCP4](https://www.ccp4.ac.uk/download/) to ensure refmac5 is accessible to the program. 

##### 4. Install LocScale and dependencies using pip:

###### Recommended installation
We recommend using pip for installation. Use pip version 21.3 or later to ensure all packages and their version requirements are met. 

```bash
pip install locscale 
```

###### Install development version
If you would like to install the latest development version of locscale, use the following command to install from the git repository. 
```bash
pip install git+https://gitlab.tudelft.nl/aj-lab/locscale.git
```

To install the git repository in editable mode, clone the repository, navigate to the `locscale` directory, and run `pip install -e .`

##### 5. Testing

To test functionality after installation, you can run LocScale unit tests using the following command:

```bash
locscale test
```

## How to use

LocScale can generate locally sharpened cryo-EM maps either using model-based sharpening based on available atomic model(s), using model-free sharpening, or using confidence-aware deep neural network-based density modification method (EMmerNet).

#### 1. Run LocScale using an existing atomic model:

```bash
locscale -em path/to/emmap.mrc -mc path/to/model.pdb -res 3 -v -o model_based_locscale.mrc
```

Here, emmap.mrc should be the unsharpened and unfiltered density map. If you wish to use the two half maps instead, use the following command:

```bash
locscale -hm path/to/halfmap1.mrc path/to/halfmap2.mrc -mc path/to/model.pdb -res 3 -v -o model_based_locscale.mrc
```

The output will be a locally sharpened map scaled according to the refined atomic B-factor distribution of the supplied atomic model.

To speed up computation, you can use multiple CPUs if available. LocScale uses [OpenMPI](https://www.open-mpi.org/)/[`mpi4py`](https://mpi4py.readthedocs.io/en/stable/) for parallelisation, which should have been automatically set up during installation. You can run it as follows:

```bash
mpirun -np 4 locscale -em path/to/emmap.mrc -mc path/to/model.pdb -res 3 -v -o model_based_locscale.mrc -mpi
```
#### 2. Run LocScale using a partial atomic model:

```bash
locscale -em path/to/emmap.mrc -mc path/to/model.pdb -res 3 -v -o model_based_locscale.mrc --complete_model
```

Here, emmap.mrc should be the unsharpened and unfiltered density map. If you wish to use the two half maps instead, use the following command:

```bash
locscale -hm path/to/halfmap1.mrc path/to/halfmap2.mrc -mc path/to/model.pdb -res 3 -v -o model_based_locscale.mrc --complete_model
```
##### Symmetry
If your map has point group symmetry, you need to specify the symmetry to force the pseudomodel generator for produce a symmetrised reference map for scaling. You can do this by specifying the required point group symmetry using the `-sym/--symmetry` flag, e.g. for D2:

```bash
locscale -hm path/to/halfmap1.mrc path/to/halfmap2.mrc -mc path/to/model.pdb -res 3 -v -sym D2 -o model_based_locscale.mrc --complete_model 
```

The output will be a locally sharpened map scaled according to the refined atomic B-factor distribution of the supplied atomic model.

To speed up computation, you can use multiple CPUs if available. LocScale uses [OpenMPI](https://www.open-mpi.org/)/[`mpi4py`](https://mpi4py.readthedocs.io/en/stable/) for parallelisation, which should have been automatically set up during installation. You can run it as follows:

```bash
mpirun -np 4 locscale -em path/to/emmap.mrc -mc path/to/model.pdb -res 3 -v -o model_based_locscale.mrc  --complete_model -mpi
```

#### 3. Run LocScale without atomic model:

If no atomic model is available, or if you do not want to use prior model information, you can use the model-free mode of `LocScale`. This method will predict a reference map using the `EMmerNet` network by default. 

Another option would be to use pseudo-atomic model. This can be enabled by passing the `--build_using_pseudomodel` flag. This mode will estimate the molecular volume using statistical thresholding and generate a pseudo-atomic model in the thresholded density map to approximate the distribution of atomic scatterers and estimate the local B-factor. It will then generate an average reference profile for local sharpening based on the experimental data and expected properties for electron scattering of biological macromolecules [[2]](#references). Use this if the default EMmerNet-based reference map generation does not work well for your data (e.g. if the map is too noisy or if the map has very low resolution). 

Usually all default parameters for pseudomodel and reference profile generation are fine, but you can [change](https://gitlab.tudelft.nl/aj-lab/locscale/-/wikis/home/) them if you deem fit.

```bash
locscale -em path/to/emmap.mrc -res 3 -v -o model_free_locscale.mrc
```
##### Symmetry
If your map has point group symmetry, you need to specify the symmetry to force the pseudomodel generator for produce a symmetrised reference map for scaling. You can do this by specifying the required point group symmetry using the `-sym/--symmetry` flag, e.g. for D2:

```bash
locscale -em path/to/emmap.mrc -res 3 -sym D2 -v -o model_free_locscale.mrc
```

LocScale currently supports all common point group symmetries. We are working on supporting helical symmetry, but this is not yet implemented. 

For faster computation, use [OpenMPI](https://www.open-mpi.org/):

```bash
mpirun -np 4 locscale -em path/to/emmap.mrc -res 3 -sym D2 -v -o model_free_locscale.mrc -mpi
```


For an exhaustive list of options, use:   

```bash
locscale --help
``` 

Alternatively, see [here](https://gitlab.tudelft.nl/aj-lab/locscale/-/wikis/home/) for more information. Please note that these pages are still being updated.

#### 4. Confidence-aware density modification:

Instead of a reference-based sharpening procedure, LocScale also supports density modification based on a physics-inspired deep neural network prediction method using  `EMmerNet` that is under development. To mitigate the risk of network inpainting or hallucination, we also calculate a per-pixel confidence score which informs the user how much the predicted density estimate deviates from a pure amplitude based sharpening approach. This score can be accessed using the p_value_map.mrc output file. We strongly recommend the users to use this score to validate the results of the prediction. 

 While we encourage its use for model building and visualisation, we do not recommend using the prediction from the neural network as target for model refinement.

```bash
locscale feature_enhance -em path/to/emmap.mrc -v -gpus 1 -o feature_enhanced_prediction.mrc
```
This will output a feature enhanced map together with a p-value map that can be used to assess the quality of the prediction.

A network trained against high context data is used as default. This will result in a prediction that shows even the weakest features in the map as long as the signal is statistically better than the noise (for example, detergent belts for membrane proteins). If you want to use a network trained against low context data, use the flag `--use_low_context_model`:

Additional models may become available and will be listed here.

For an exhaustive list of options, run:   

```bash
locscale feature_enhance --help
``` 

## Tutorial and FAQs

We are currently working on the tutorial and [__Wiki__](https://gitlab.tudelft.nl/aj-lab/locscale/-/wikis/home/) help. If you are still using LocScale 1.0, see the [LocScale1](https://gitlab.tudelft.nl/ajakobi/locscale/wikis/home)-Wiki for usage instructions, FAQs and tutorial.
<br>  

## Credits

This project is using code from a number of third-party open-source projects. Projects used by `LocScale` are included under include/:

[EMmer](https://gitlab.tudelft.nl/aj-lab/emmer) - Python library for electron microscopy map and model manipulations. License: 3-Clause BSD.     
[FDRthresholding](https://git.embl.de/mbeckers/FDRthresholding) – tool for FDR-based density thresholding. License: 3-Clause BSD.     

`LocScale` also makes use of [Refmac](https://www2.mrc-lmb.cam.ac.uk/groups/murshudov/content/refmac/refmac.html) – coordinate refinement program for macromolecular structures. Refmac is distributed as part of CCP-EM.

## References

If you found `LocScale` useful for your research, please consider citing it:

- A.J. Jakobi, M. Wilmanns and C. Sachse, [Model-based local density sharpening of cryo-EM maps](https://doi.org/10.7554/eLife.27131), eLife 6: e27131 (2017).
- A. Bharadwaj and A.J. Jakobi, [Electron scattering properties and their use in cryo-EM map sharpening](https://doi.org/10.1039/D2FD00078D), Faraday Discussions 240, 168-183 (2022)
---

## Bugs and questions

For bug reports please use the [GitLab issue tracker](https://gitlab.tudelft.nl/aj-lab/locscale/issues).   
