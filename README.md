# landscape_physiographic_patchness
This repository stores MATLAB codes of a patch model that calculates and visualizes physiographic patches of a landscape evolution model and patch changes through time.  

This habitat model was invented to explore the habitat (combinations of physiographic parameters' value ranges) evolution of a retreating escarpment. The manuscript _"Escarpment Evolution Drives the Diversification of the Madagascar Flora"_ has been submitted to _Science_. Revisions of the manuscript will be submitted soon. 

## 1. Physiographic patch
The land surface (location and height) can be represented by a series of elemental surfaces ( _Si(x, y, z)_ ). When defining a certain range of physiographic parameters, for example, elevation or hillslope gradient, elemental surfaces with parametric values falling in this range may be connected and form a contiguous patch of surface. Physiographic parameters used in this model are elevation, hillslope gradient, and hillslope aspect.

The figure below shows a modeled landscape (colored), with patches (dark) that have elevation between 500-1500 m, hillslope gradient of 3-20 degrees, and aspect not facing north.  
![escarpment_topo_example_patch_3](https://github.com/yanyanwangesd/landscape_physiographic_patchness/assets/108676831/2aaf2b8c-15e6-4608-a60e-332218d99619)

## 2. Patch model inputs
The patch model uses topographic data (location and height, _x_, _y_, _z_ ) of landscape evolution models (LEMs) to calculate patches that are defined by physiographic parameters. 

## 3. Numerical definition of patches
Topographic data from LEMs are read into the MATLAB code. The topography data are interpolated onto a regular grid when necessary because some LEMs use irregular grids. Physiographic parameters are calculated for elemental surfaces. When defining a combination of ranges of physiographic parameters, elemental surfaces that meet this combination would form clusters on the (_x_, _y_) plane, shown as the dark cluster on the figure below. These patches are then identified with image functions in MATLAB.  
![escarpment_topo_example_patch](https://github.com/yanyanwangesd/landscape_physiographic_patchness/assets/108676831/8f08157a-02d8-4e97-a52b-2a966415c59d)

## 4. Patch change types
Each individual patch is tracked and flagged with unique identities to track the temporal changes of patches. Types of patch change are birth, death, split, merge, and deform. Consecutive steps are compared to identify different types of patch changes. 
![patch_evolution_carton_illustration_version2](https://github.com/yanyanwangesd/landscape_physiographic_patchness/assets/108676831/32f4828c-cbf5-4f8e-857f-258fd84af6cf)




