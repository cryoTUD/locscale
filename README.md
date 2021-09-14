[![PyPI - Python Version](https://img.shields.io/pypi/pyversions/locscale)](https://pypi.org/project/locscale)
[![PyPI](https://img.shields.io/pypi/v/instamatic.svg?style=flat)](https://pypi.org/project/locscale/)
[![PyPI - Downloads](https://img.shields.io/pypi/dm/locscale)](https://pypi.org/project/locscale/)
[![DOI](https://zenodo.org/badge/DOI/10.7554/eLife.2713110.1007.svg)](https://doi.org/10.7554/eLife.27131)

# LocScale

`LocScale` is a program for reference-based density scaling (local sharpening) of cryo-EM maps with the aim to improve their interpretability.
  
`LocScale` is distributed as a portable stand-alone installation that includes all the needed libraries from: https://gitlab.tudelft.nl/aj-lab/locscale/releases.  

Please note that there is a GUI implemented version available as part of the [CCP-EM](http://www.ccpem.ac.uk/) project; it is also implemented in [Scipion](http://scipion.i2pc.es/).
<br>   

## Installation

If you use conda, create a new environment:

```
conda create -n locscale python=3.8
conda activate locscale
```

Install using pip:

```bash
pip install locscale
```

Alternatively, download the portable installation with all libraries/dependencies included: https://gitlab.tudelft.nl/aj-lab/locscale/releases/latest.
<br> 

## Usage, tutorial and FAQs

Please see the [__Wiki__](https://gitlab.tudelft.nl/ajakobi/locscale/wikis/home) pages for usage instructions, FAQs and tutorial.
<br>  

## Credits

This project is using code from a number of third-party open-source projects. Projects used in the C++ library and included under include/:

[FDRthresholding](https://git.embl.de/mbeckers/FDRthresholding) – tool for FDR-based density thresholding. License: BSD 3-clause.  
[EMDA](https://www2.mrc-lmb.cam.ac.uk/groups/murshudov/content/emda/emda.html) – Python library for electron microscopy map and model manipulations. Licence: MPL-2

`LocScale` also makes use of [Refmac](https://www2.mrc-lmb.cam.ac.uk/groups/murshudov/content/refmac/refmac.html) – coordinate refinement program for macromolecular structures. Refmac is distributed as part of CCP-EM.

## Reference

If you found `LocScale` useful, please consider citing it:

- A. Bharadwaj and A.J. Jakobi, [Title still to be determined](do-not-know-yet), A Journal 1, 1-2 (2021)
- A.J. Jakobi, M. Wilmanns and C. Sachse, [Model-based local density sharpening of cryo-EM maps](https://doi.org/10.7554/eLife.27131), eLife 6: e27131 (2017).

---


For bug reports please use the [GitLab issue tracker](https://gitlab.tudelft.nl/aj-lab/locscale/issues).   
For questions or comments please contact <a.jakobi@tudelft.nl> or a.bharadwaj@tudelft.nl.
