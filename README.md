# landscape_physiographic_patchness
[![DOI](https://zenodo.org/badge/684976091.svg)](https://zenodo.org/badge/latestdoi/684976091)

This repository stores MATLAB codes of a patch model that calculates and visualizes physiographic patches of a landscape evolution model and patch changes through time.  

This landscape physiographic patch model was invented to explore the characteristic habitat (combinations of physiographic parameters' value ranges) evolution of a retreating escarpment. The manuscript _"Escarpment Evolution Drives the Diversification of the Madagascar Flora"_ has been submitted to _Science_. Revisions of the manuscript will be submitted soon. 

## 1. Physiographic patch
The land surface (location and height) can be represented by a series of elemental surfaces ( _Si(x, y, z)_ ). When defining a certain range of physiographic parameters, for example, elevation or hillslope gradient, elemental surfaces with parametric values falling in this range may be connected and form a contiguous patch of surface. Physiographic parameters used in this model are elevation, hillslope gradient, and hillslope aspect.

The figure below shows a modeled landscape (colored), with patches (dark) that have elevation between 500-1500 m, hillslope gradient of 3-20 degrees, and aspect not facing north.  
![escarpment_topo_example_patch_3](https://github.com/yanyanwangesd/landscape_physiographic_patchness/assets/108676831/2aaf2b8c-15e6-4608-a60e-332218d99619)

## 2. Patch model inputs
The patch model uses topographic data (location and height, _x_, _y_, _z_ ) of landscape evolution models (LEMs) to calculate patches that are defined by physiographic parameters. 

## 3. Numerical definition of patches
Topographic data from LEMs are read into the MATLAB code. The topography data are interpolated onto a regular grid when necessary because some LEMs use irregular grids. Physiographic parameters are calculated for elemental surfaces. When defining a combination of ranges of physiographic parameters, elemental surfaces that meet this combination would form clusters on the (_x_, _y_) plane, shown as the dark cluster in the figure below. These patches are then identified with image functions in MATLAB.  
![escarpment_topo_example_patch](https://github.com/yanyanwangesd/landscape_physiographic_patchness/assets/108676831/8f08157a-02d8-4e97-a52b-2a966415c59d)

## 4. Patch change types
Each individual patch is tracked and flagged with unique identities to track the temporal changes of patches. Types of patch change are birth, death, split, merge, and deform. Consecutive steps are compared to identify different types of patch changes. 
![patch_evolution_carton_illustration_version2](https://github.com/yanyanwangesd/landscape_physiographic_patchness/assets/108676831/32f4828c-cbf5-4f8e-857f-258fd84af6cf)


## 5. What's in the repository?
- The ASCII folder stores example input files.
- The script _landscape_physiographic_patchness.m_ calculates patches.
- The script _plot_patch_change_histograms.m_ plots histograms of patch changes. 

## 6. Technical help
Technical support or discussions of the methodology are welcome, please contact Dr. Yanyan Wang via email (yanyan.wang@erdw.ethz.ch or wangyanyan0607@hotmail.com).

## 7. Citation
- Please cite the paper Liu, Y., Wang, Y., Willett, S.D., Zimmermann, N.E. and Pellissier, L., 2024. Escarpment evolution drives the diversification of the Madagascar flora. Science, 383(6683), pp.653-658. DOI: [10.1126/science.adi0833](https://doi.org/10.1126/science.adi0833)
- Please cite the repository through [![DOI](https://zenodo.org/badge/684976091.svg)](https://zenodo.org/badge/latestdoi/684976091)
- The paper is free-access through the URL: [https://www.science.org/stoken/author-tokens/ST-1688/full](https://www.science.org/stoken/author-tokens/ST-1688/full)

  Bibtex:
```
@article{liu2024escarpment,
  title={Escarpment evolution drives the diversification of the Madagascar flora},
  author={Liu, Yi and Wang, Yanyan and Willett, Sean D and Zimmermann, Niklaus E and Pellissier, Lo{\"\i}c},
  journal={Science},
  volume={383},
  number={6683},
  pages={653--658},
  year={2024},
  publisher={American Association for the Advancement of Science}
}
```


