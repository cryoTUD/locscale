[![PyPI - Python Version](https://img.shields.io/pypi/pyversions/locscale)](https://pypi.org/project/locscale)
[![PyPI](https://img.shields.io/pypi/v/locscale.svg?style=flat)](https://pypi.org/project/locscale/)
[![PyPI - Downloads](https://img.shields.io/pypi/dm/locscale)](https://pypi.org/project/locscale/)
[![DOI](https://zenodo.org/badge/DOI/10.7554/eLife.2713110.1007.svg)](https://doi.org/10.7554/eLife.27131)

# LocScale - reference-based density scaling for local sharpening of cryo-EM maps

`LocScale` is an automated program for local sharpening of cryo-EM maps with the aim to improve their interpretability. It utilises general properties inherent to electron scattering from biological macromolecules to restrain the sharpening filter. These can be provided either from an existing atomic model, or inferred directly from the experimental density map.

#### New in LocScale 2.0:

- Model-free sharpening: `LocScale` now supports reference-based sharpening without the need to supply an atomic model

- Completely automated process for local map sharpening 

- Full support for point group symmetry (helical symmetry to follow)

- `EMmerNet`: deep convolutional neural network-based sharpening method. `EMmerNet` is an ensemble network model trained on model-free `LocScale` maps from a large number of existing cryo-EM structures in the [EMDB]().
<br>
  
`LocScale` is distributed as a portable stand-alone installation that includes all the needed libraries from: https://gitlab.tudelft.nl/aj-lab/locscale/releases.   


Please note that there is a GUI implemented version available as part of the [CCP-EM](http://www.ccpem.ac.uk/) project; it is also implemented in [Scipion](http://scipion.i2pc.es/).
<br>   

## Installation 

We recommend to use [Conda](https://docs.conda.io/en/latest/) for a local working environment. See [here](https://docs.conda.io/projects/conda/en/latest/user-guide/install/download.html#anaconda-or-miniconda) for more information on what Conda flavour may be the right choice for you, and [here](https://www.anaconda.com/products/distribution) for Conda installation instructions.

#### 1. Create a and acticate a new conda environment

```bash
conda create -n locscale python=3.8 
source activate locscale
conda install -c conda-forge gfortran
```

#### 2. Install LocScale and dependencies using pip:

The setup.py file contains the list of packages and their versions used inside LocScale. Use pip version 21.3 or later to ensure all packages and their version requirements are met. 

```bash
pip install locscale 
```

Alternatively, download the portable installation with all libraries/dependencies included: https://gitlab.tudelft.nl/aj-lab/locscale/releases/latest.


## Usage

#### 1. Run LocScale using an existing atomic model:
```bash
locscale run_locscale -em path/to/emmap.mrc -mc path/to/model.pdb -res 3 -v -o model_based_locscale.mrc
```

#### 2. Run LocScale without atomic model:
```bash
locscale run_locscale -em path/to/emmap.mrc -res 3 -v -o model_based_locscale.mrc
```

For an exhaustive list of options, run:   

```bash
locscale run_locscale --help
``` 

#### 3. Run LocScale using EMmerNet predictions:
```bash
locscale run_emmernet -em path/to/emmap.mrc -v -trained_model model_based -gpus 0 -o emmernet_model_based.mrc
```

Currently, three different EMmerNet models are available and can be specified using the `-trained_model` flag as follows:

| Model  | Syntax  | 
|---|---|
| Model Based:       | ```-trained_model model_based```| 
| Model Free:        | ```-trained_model model_free``` | 
| Ensemble Network:  | ```-trained_model ensemble```   | 

Additional models may become available and will be listed here.

For an exhaustive list of options, run:   

```bash
locscale run_emmernet --help
``` 


## Tutorial and FAQs

We are currently working on the tutorial and Wiki help. If you are still using LocScale 1.0, see [__Wiki__](https://gitlab.tudelft.nl/ajakobi/locscale/wikis/home) for usage instructions, FAQs and tutorial.
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
For questions or comments please contact <a.jakobi@tudelft.nl> or a.bharadwaj@tudelft.nl.
