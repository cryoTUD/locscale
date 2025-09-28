<!---[![stability-beta](https://img.shields.io/badge/stability-beta-33bbff.svg)](https://github.com/mkenney/software-guides/blob/master/STABILITY-BADGES.md#beta)-->
[![stability-release-candidate](https://img.shields.io/badge/stability-pre--release-48c9b0.svg)](https://github.com/mkenney/software-guides/blob/master/STABILITY-BADGES.md#release-candidate)
[![Python 3.11](https://img.shields.io/badge/python-3.11-green)](https://www.python.org/downloads/release/python-3110/)
[![PyPI](https://img.shields.io/pypi/v/locscale.svg?style=flat)](https://pypi.org/project/locscale/)
[![PyPI - Downloads](https://img.shields.io/pypi/dm/locscale)](https://pypi.org/project/locscale/)
[![License](https://img.shields.io/pypi/l/locscale.svg?color=orange)](https://gitlab.tudelft.nl/aj-lab/locscale/raw/master/LICENSE)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.15488220.svg)](https://doi.org/10.5281/zenodo.15488220)
[![Citation Badge](https://api.juleskreuer.eu/citation-badge.php?doi=10.7554/eLife.27131)](https://juleskreuer.eu/projekte/citation-badge/)

# LocScale-2.0
`LocScale-2.0` is an automated map optimisation program performing physics-informed local sharpening and/or density modification of cryo-EM maps with the aim to improve their interpretability. It utilises general properties inherent to electron scattering from biological macromolecules to restrain the sharpening and/or optimisation filter. These can be inferred directly from the experimental density map, or - in legacy mode – provided from an existing atomic model.

## What's new in LocScale-2.0?
- Completely automated process for map optimisation
- [`Feature-enhanced maps`](https://locscale.readthedocs.io/en/latest/tutorials/fem/): Confidence-weighted map optimisation by variational inference.
- [`Hybrid sharpening`](https://locscale.readthedocs.io/en/latest/tutorials/hybrid_locscale/): Reference-based local sharpening with partial (incomplete) models. 
- [`Model-free sharpening`](https://locscale.readthedocs.io/en/latest/tutorials/model_free_locscale/): Reference-based local sharpening without atomic models.
- [`LocScale-SURFER`](https://github.com/cryoTUD/locscale-surfer): ChimeraX plugin to toggle contextual structure in LocScale maps.
- Full support for point group symmetry (helical symmetry to follow).

## Documentation

>[!IMPORTANT]
> Please visit [https://cryotud.github.io/locscale/](https://cryotud.github.io/locscale/) for comprehensive documentation, tutorials and troubleshooting.

## Installation

We recommend to use [Conda](https://docs.conda.io/en/latest/) for a local working environment. See [here](https://docs.conda.io/projects/conda/en/latest/user-guide/install/download.html#anaconda-or-miniconda) for more information on what Conda flavour may be the right choice for you, and [here](https://www.anaconda.com/products/distribution) for Conda installation instructions.

>[!NOTE]
>LocScale should run on any CPU system with Linux, OS X or Windows subsytem for Linux (WSL). To run LocScale efficiently in EMmerNet mode requires the availability of a GPU; it is possible to run it on CPUs but computation will be slow(er). 

<details>

<summary> Quick installation </summary>

We recommend to use [Conda](https://docs.conda.io/en/latest/) for a local working environment. See [here](https://docs.conda.io/projects/conda/en/latest/user-guide/install/download.html#anaconda-or-miniconda) for more information on what Conda flavour may be the right choice for you, and [here](https://www.anaconda.com/products/distribution) for Conda installation instructions.

##### 1. Install `LocScale-2.0` using environment files 

Download [environment.yml](https://github.com/cryoTUD/locscale/blob/master/environment.yml) to your local computer, navigate to the location you wish to install `Locscale-2.0` at and run the following: 

```bash
conda env create -f /path/to/environment.yml
conda activate locscale2
```

##### 2. Install REFMAC5 via CCP4/CCPEM
`LocScale` needs a working instance of [REFMAC5](https://www2.mrc-lmb.cam.ac.uk/groups/murshudov/index.html). If you already have CCP4/CCPEM installed check if the path to run `refmac5` is present in your environment. 

```bash
which refmac5
```

If no valid path is returned, please install [CCP4](https://www.ccp4.ac.uk/download/) to ensure refmac5 is accessible to the program. 
</details>

<details>
<summary> Step-by-step instructions </summary>

##### 1. Create and activate a new conda environment

```bash
conda create -n locscale python=3.11
conda activate locscale
```

##### 2. Install parallelisation support and Fortran compiler
`LocScale` uses Fortran code to perform symmetry operations and requires a Fortran compiler to be present in your system.  You can install `gfortran`, `mpi4py` and `openmpi` from conda-forge. 
```bash
conda install -c conda-forge gfortran mpi4py openmpi
```
##### 3. Install REFMAC5 via CCP4/CCPEM

The model-based and hybrid map sharpening modes of LocScale need a working instance of [REFMAC5](https://www2.mrc-lmb.cam.ac.uk/groups/murshudov/index.html). If you already have CCP4/CCPEM installed check if the path to run `refmac5` is present in your environment. For model-free sharpening and confidence-aware density modification REFMAC5 is not required. 

```bash
which refmac5
```

If no valid path is returned, please install [CCP4](https://www.ccp4.ac.uk/download/) to ensure REFMAC5 is accessible to the program. 

##### 4. Install LocScale and dependencies using pip:

We recommend using pip for installation. Use pip version 21.3 or later to ensure all packages and their version requirements are met. 

```bash
pip install locscale 
```

>[!NOTE] 
> ##### Install development version: 
>If you would like to install the latest development version of locscale, use the following command to install from the git repository. 
>```bash
>pip install git+https://github.com/cryoTUD/locscale.git
>```

To install the git repository in editable mode, clone the repository, navigate to the `locscale` directory, and run `pip install -e .`

##### 5. Testing

To test functionality after installation, you can run LocScale unit tests using the following command:

```bash
locscale test
```
</details>

## ColabScale

>[!TIP]
> For quick testing or if you have limited compute resources, many functionalities of `LocScale-2.0` are available on [`ColabScale`](https://github.com/cryoTUD/ColabScale/tree/main).  <br><br>
> <a target="_blank" href="https://colab.research.google.com/github/cryoTUD/ColabScale/blob/main/ColabScale.ipynb">
>  <img src="https://colab.research.google.com/assets/colab-badge.svg" alt="Open In Colab"/>
>  </a>

## Credits
`LoScale 2.0` is facilitated by a number of open-source projects.

- [`EMmer`](https://gitlab.tudelft.nl/aj-lab/emmer): Python library for electron microscopy map and model manipulations. [3-Clause BSD license]    
- [`FDRthresholding`](https://git.embl.de/mbeckers/FDRthresholding): Tool for FDR-based density thresholding. [3-Clause BSD license]
- [`EMDA`](https://gitlab.com/ccpem/emda/): Electron Microscopy Data Analytical Toolkit. [MPL2.0 license]
- [`Servalcat`](https://github.com/keitaroyam/servalcat): Structure refinement and validation for crystallography and SPA. [MPL2.0 license]
- [`mrcfile`](https://pypi.org/project/mrcfile/): MRC file I/O. [3-Clause BSD license]

`LocScale` also makes use of [REFMAC5](https://www2.mrc-lmb.cam.ac.uk/groups/murshudov/content/refmac/refmac.html). REFMAC is distributed as part of CCP-EM.

## References

If you found `LocScale` useful for your research, please consider citing it:

- A. Bharadwaj, R.M. de Bruin, A.J. Jakobi: [Confidence-guided cryo-EM map optimisation with LocScale-2.0](https://doi.org/10.1101/2025.09.11.674726), BioRxiv 2025.09.11.674726 (2025) 
- A.J. Jakobi, M. Wilmanns and C. Sachse: [Model-based local density sharpening of cryo-EM maps](https://doi.org/10.7554/eLife.27131), eLife 6: e27131 (2017).
- A. Bharadwaj and A.J. Jakobi: [Electron scattering properties and their use in cryo-EM map sharpening](https://doi.org/10.1039/D2FD00078D), Faraday Discussions 240, 168-183 (2022)
---

## Bugs and questions

For bug reports please use the [GitHub issue tracker](https://github.com/issues/assigned).   
