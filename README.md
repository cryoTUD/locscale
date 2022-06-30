# Development
This branch is currently under development. Please switch to the master branch for a working installation. 

# LocScale using predicted profiles
The goal of this project is to perform Local sharpening using predicted profiles. 

# Authors:
- Ian Bot, (TU Delft)
- Alok Bharadwaj, (TU Delft)
- Arjen Jakobi, (TU Delft)

[![Python 3.6](https://img.shields.io/badge/python-3.7%20%7C%203.8-brightgreen)](https://www.python.org/downloads/release/python-370/)
[![PyPI](https://img.shields.io/pypi/v/locscale.svg?style=flat)](https://pypi.org/project/locscale/)
[![PyPI - Downloads](https://img.shields.io/pypi/dm/locscale)](https://pypi.org/project/locscale/)
[![License](https://img.shields.io/pypi/l/locscale.svg?color=orange)](https://gitlab.tudelft.nl/aj-lab/locscale/raw/master/LICENSE)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.6652013.svg)](https://doi.org/10.5281/zenodo.6652013)
[![Citation Badge](https://api.juleskreuer.eu/citation-badge.php?doi=10.7554/eLife.27131)](https://juleskreuer.eu/projekte/citation-badge/)

# LocScale - reference-based local sharpening of cryo-EM maps

`LocScale` is an automated program for local sharpening of cryo-EM maps with the aim to improve their interpretability. It utilises general properties inherent to electron scattering from biological macromolecules to restrain the sharpening filter. These can be provided either from an existing atomic model, or inferred directly from the experimental density map.

## Installation 

#### Requirements

LocScale should run on any CPU system with Linux, OS X or Windows subsytem for Linux (WSL). To run LocScale efficiently in EMmerNet mode requires the availability of a GPU; it is possible to run it on CPUs but computation will be slow. 

#### Installation instructions:

##### 1. Create and activate a new conda environment

```bash
conda create -n locscale python=3.8 
conda activate locscale
```

###### Install development version
If you would like to install the latest development version of locscale, use the following command to install from the git repository. 
```bash
pip install git+https://gitlab.tudelft.nl/aj-lab/locscale.git@predict_profile
```

To install the git repository in editable mode, clone the repository and navigate to the `locscale` directory. Switch branch:  `git checkout predict_profile` and run `pip install -e .`

## References

If you found `LocScale` useful, please consider citing it:

- A.J. Jakobi, M. Wilmanns and C. Sachse, [Model-based local density sharpening of cryo-EM maps](https://doi.org/10.7554/eLife.27131), eLife 6: e27131 (2017).
- A. Bharadwaj and A.J. Jakobi, [Electron scattering properties and their use in cryo-EM map sharpening](https://doi.org/10.1039/D2FD00078D), Faraday Discussions D2FD00078D (2022)
---

## Bugs and questions

For bug reports please use the [GitLab issue tracker](https://gitlab.tudelft.nl/aj-lab/locscale/issues).   
