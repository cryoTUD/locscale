`LocScale 2.0` has been developed by Alok Bharadwaj and Arjen Jakobi at the [Kavli Institute of Nanoscience Delft](https://kavli.tudelft.nl/), with contributions by Reinier de Bruijn to the development of EMmerNet. 
<br>

### Credits
`LoScale 2.0` is facilitated by a number of open-source projects.

!!! info "Projects used in LocScale"

    - [`EMmer`](https://gitlab.tudelft.nl/aj-lab/emmer): Python library for electron microscopy map and model manipulations. [3-Clause BSD license]    
    - [`FDRthresholding`](https://git.embl.de/mbeckers/FDRthresholding): Tool for FDR-based density thresholding. [3-Clause BSD license]
    - [`EMDA`](https://gitlab.com/ccpem/emda/): Electron Microscopy Data Analytical Toolkit. [MPL2.0 license]
    - [`Servalcat`](https://github.com/keitaroyam/servalcat): Structure refinement and validation for crystallography and SPA. [MPL2.0 license]
    - [`mrcfile`](https://pypi.org/project/mrcfile/): MRC file I/O. [3-Clause BSD license]

    `LocScale` also makes use of [REFMAC5](https://www2.mrc-lmb.cam.ac.uk/groups/murshudov/content/refmac/refmac.html). Refmac is distributed as part of CCP-EM.

We acknowledge the use of computational resources of the DelftBlue supercomputer, provided by [Delft High Performance Computing Centre](https://www.tudelft.nl/dhpc).  

## References {#references}

If you found `LocScale-2.0` useful for your research, please consider citing it:

If you found `LocScale` useful for your research, please consider citing it:
```
@article{Bharadwaj2026,
  author  = {Bharadwaj, A. and de Bruin, R. M. and Jakobi, A. J.},
  title   = {Confidence-guided cryo-EM map optimisation with LocScale-2.0},
  journal = {Nature Communications},
  volume  = {17},
  pages   = {8778},
  year    = {2026},
  doi     = {10.1038/s41467-026-75327-8}
}
```
```
@article{Jakobi2017,
  author  = {Jakobi, A. J. and Wilmanns, M. and Sachse, C.},
  title   = {Model-based local density sharpening of cryo-EM maps},
  journal = {eLife},
  volume  = {6},
  pages   = {e27131},
  year    = {2017},
  doi     = {10.7554/eLife.27131}
}
```


---

## Bugs and questions

For bug reports please use the [LocScale issue tracker](https://github.com/cryoTUD/locscale/issues). 
