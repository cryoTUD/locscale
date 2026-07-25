<!---[![stability-beta](https://img.shields.io/badge/stability-beta-33bbff.svg)](https://github.com/mkenney/software-guides/blob/master/STABILITY-BADGES.md#beta)-->
[![stability-release-candidate](https://img.shields.io/badge/stability-pre--release-48c9b0.svg)](https://github.com/mkenney/software-guides/blob/master/STABILITY-BADGES.md#release-candidate)
[![Python 3.12](https://img.shields.io/badge/python-3.12-green)](https://www.python.org/downloads/release/python-3120/)
[![PyPI](https://img.shields.io/pypi/v/locscale.svg?style=flat)](https://pypi.org/project/locscale/)
[![PyPI - Downloads](https://img.shields.io/pypi/dm/locscale)](https://pypi.org/project/locscale/)
[![License](https://img.shields.io/pypi/l/locscale.svg?color=orange)](https://github.com/cryoTUD/locscale/blob/master/LICENSE)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.15488220.svg)](https://doi.org/10.5281/zenodo.15488220)
[![Citations LocScale](https://api.juleskreuer.eu/citation-badge.php?doi=10.7554/eLife.27131)](https://doi.org/10.7554/eLife.27131)
[![Citations LocScale2](https://api.juleskreuer.eu/citation-badge.php?doi=10.1038/s41467-026-75327-8)](https://www.nature.com/articles/s41467-026-75327-8)

# LocScale Feature Enhance for ChimeraX


## What's new?
- You can run LocScale Feature Enhance directly in your ChimeraX workspace! 

## Documentation

>[!IMPORTANT]
> Please visit [https://cryotud.github.io/locscale/](https://cryotud.github.io/locscale/) for documentation about LocScale. 

## Installation
- Clone this branch
```bash
git clone https://github.com/cryoTUD/locscale.git@locscale_fem_chimerax
```

- Inside your ChimeraX command line: 
```chimerax
devel install path\to\locscale\
```

or (coming soon!)
```chimerax
toolshed install locscale2
```

## Usage

- Verify LocScale2 figures from the paper. This is useful for a quick check of the published results

```chimerax
locscale2 verify 33888
```
It downloads the ChimeraX session containing the maps and views present in the paper. It also downloads the halfmaps from the EMDB, symmetry information and preloads everything into the LocScale2 GUI on ChimeraX. Click run feature enhance to check the results. On the log you should see the correlation of the predicted map with the published map. 

```chimerax
locscale2 verify list
```

shows all EMDBs that can be verified through published chimerax sessions. Currently EMDB IDs in Figure 7/Supplementary 8 can be verified through this tool. 


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
