[![PyPI - Python Version](https://img.shields.io/pypi/pyversions/locscale)](https://pypi.org/project/locscale)
[![PyPI](https://img.shields.io/pypi/v/locscale.svg?style=flat)](https://pypi.org/project/locscale/)
[![PyPI - Downloads](https://img.shields.io/pypi/dm/locscale)](https://pypi.org/project/locscale/)
[![DOI](https://zenodo.org/badge/DOI/10.7554/eLife.2713110.1007.svg)](https://doi.org/10.7554/eLife.27131)

# LocScale - reference-based local sharpening of cryo-EM maps

`LocScale` is an automated program for local sharpening of cryo-EM maps with the aim to improve their interpretability. It utilises general properties inherent to electron scattering from biological macromolecules to restrain the sharpening filter. These can be provided either from an existing atomic model, or inferred directly from the experimental density map.

#### New in LocScale 2.0:

- Model-free sharpening: `LocScale` now supports reference-based sharpening without the need to supply an atomic model

- Completely automated process for local map sharpening 

- Full support for point group symmetry (helical symmetry to follow)

- `EMmerNet`: deep convolutional neural network-based sharpening method. `EMmerNet` is an ensemble network model trained on model-free LocScale maps from a large number of existing cryo-EM structures in the [EMDB]().
<br>
  
`LocScale` is distributed as a portable stand-alone installation that includes all the needed libraries from: https://gitlab.tudelft.nl/aj-lab/locscale/releases.   


Please note that there is a GUI implemented version available as part of the [CCP-EM](http://www.ccpem.ac.uk/) project; it is also implemented in [Scipion](http://scipion.i2pc.es/). Note that currently the GUI implementations only support an older version of Locscale (Locscale 1.0, with only model-based sharpening). 

## Installation 

We recommend to use [Conda](https://docs.conda.io/en/latest/) for a local working environment. See [here](https://docs.conda.io/projects/conda/en/latest/user-guide/install/download.html#anaconda-or-miniconda) for more information on what Conda flavour may be the right choice for you, and [here](https://www.anaconda.com/products/distribution) for Conda installation instructions.

#### Requirements

LocScale should run on any CPU system with Linux, OS X or Windows subsytem for Linux (WSL). To run LocScale efficiently in EMmerNet mode requires the availability of a GPU; it is possible to run it on CPUs but computation will be slow. 

#### Installation instructions:

##### 1. Create and activate a new conda environment

```bash
conda create -n locscale python=3.7 
conda activate locscale
```

###### Install fortran compiler
LocScale uses Fortran code to perform symmetry operations and requires a Fortran compiler to be present in your system. You can install `gfortran` from conda-forge.
```bash
conda install -c conda-forge gfortran
```

##### 2. Install LocScale and dependencies using pip:

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

###### Testing

To test functionality after installation, you can run LocScale unit tests using the following command:

```bash
locscale test
```


##### 3. Install REFMAC5 via CCP4/CCPEM

LocScale needs a working instance of [REFMAC5](https://www2.mrc-lmb.cam.ac.uk/groups/murshudov/index.html). If you already have CCP4/CCPEM installed, you can skip this step. Check if the path to run `refmac5` is present in your environment. 

```bash
which refmac5
```

If no valid path is returned, please install [CCP4](https://www.ccp4.ac.uk/download/). Note that you can provide an alternative path to REFMAC5 during runtime by using the `--refmac5_path REFMAC5_PATH` flag in LocScale.

## How to use

LocScale can generate locally sharpened cryo-EM maps either using model-based sharpening based on available atomic model(s), using model-free sharpening, or using deep neural network-based sharpening method (EMmerNet).

#### 1. Run LocScale using an existing atomic model:

```bash
locscale run_locscale -em path/to/emmap.mrc -mc path/to/model.pdb -res 3 -v -o model_based_locscale.mrc
```

Here, emmap.mrc should be the unsharpened and unfiltered density map. If you wish to use the two half maps instead, use the following command:

```bash
locscale run_locscale -hm path/to/halfmap1.mrc,path/to/halfmap2.mrc -mc path/to/model.pdb -res 3 -v -o model_based_locscale.mrc
```

The output will be a locally sharpened map scaled according to the refined atomic B-factor distribution of the supplied atomic model.

To speed up computation, you can use multiple CPUs if available. LocScale uses [OpenMPI](https://www.open-mpi.org/)/[`mpi4py`](https://mpi4py.readthedocs.io/en/stable/) for parallelisation, which should have been automatically set up during installation. You can run it as follows:

```bash
mpirun -np 4 locscale run_locscale -em path/to/emmap.mrc -mc path/to/model.pdb -res 3 -v -o model_based_locscale.mrc -mpi
```

#### 2. Run LocScale without atomic model:

If no atomic model is available, or if you do not want to use prior model information, you can use the model-free mode of `LocScale`. This mode will estimate the molecular volume using statistical thresholding and generate a pseudo-atomic model in the thresholded density map to approximate the distribution of atomic scatterers and estimate the local B-factor. It will then generate an average reference profile for local sharpening based on the experimental data and expected properties for electron scattering of biological macromolecules [[2]](#references). Usually all default parameters for pseudomodel and reference profile generation are fine, but you can [change](https://gitlab.tudelft.nl/aj-lab/locscale/-/wikis/home/) them if you deem fit.

```bash
locscale run_locscale -em path/to/emmap.mrc -res 3 -v -o model_free_locscale.mrc
```
##### Symmetry
If your map has point group symmetry, you need to specify the symmetry to force the pseudomodel generator for produce a symmetrised reference map for scaling. You can do this by specifying the required point group symmetry using the `-sym/--symmetry` flag, e.g. for D2:

```bash
locscale run_locscale -em path/to/emmap.mrc -res 3 -sym D2 -v -o model_free_locscale.mrc
```

LocScale currently supports all common point group symmetries. We are working on supporting helical symmetry, but this is not yet implemented. 

For faster computation, use [OpenMPI](https://www.open-mpi.org/):

```bash
mpirun -np 4 locscale run_locscale -em path/to/emmap.mrc -res 3 -sym D2 -v -o model_free_locscale.mrc -mpi
```


For an exhaustive list of options, use:   

```bash
locscale run_locscale --help
``` 

Alternatively, see [here](https://gitlab.tudelft.nl/aj-lab/locscale/-/wikis/home/) for more information. Please note that these pages are still being updated.

#### 3. Run LocScale using EMmerNet predictions:

Instead of using model-based or model-free refernce profiles, LocScale also supports local sharpening based on a new deep neural network prediction method using our ensemble network `EMmerNet` that is under development. While we have done our very best to validate the network and mitigate the risk of hallucination, as inherent to any such approach (including existing software) care needs to be exercised to avoid the slightest risk of bias. We do encourage its use, but ask you to be cautious and, in particular, not to use such maps for model refinement.

```bash
locscale run_emmernet -em path/to/emmap.mrc -v -trained_model model_based -gpus 1 -o emmernet_model_based.mrc
```

Currently, three different EMmerNet models are available and can be specified using the `-trained_model` flag as follows:

| Model  | Syntax  | 
|---|---|
| Model Based EMmerNet  | ```-trained_model model_based```| 
| Model Free EMmerNet   | ```-trained_model model_free``` | 
| Ensemble EMmerNet     | ```-trained_model ensemble```   | 

Additional models may become available and will be listed here.

For an exhaustive list of options, run:   

```bash
locscale run_emmernet --help
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

If you found `LocScale` useful, please consider citing it:

- A.J. Jakobi, M. Wilmanns and C. Sachse, [Model-based local density sharpening of cryo-EM maps](https://doi.org/10.7554/eLife.27131), eLife 6: e27131 (2017).
- A. Bharadwaj and A.J. Jakobi, [Electron scattering properties and their use in cryo-EM map sharpening](https://doi.org/10.1039/D2FD00078D), Faraday Discussions D2FD00078D (2022)
---

## Bugs and questions

For bug reports please use the [GitLab issue tracker](https://gitlab.tudelft.nl/aj-lab/locscale/issues).   
