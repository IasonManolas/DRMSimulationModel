This repository provides a cross-platform computational tool for predicting the static equilibrium of general bending-active structures in
the form of a network of elements using the dynamic relaxation method. The implemented method is presented in an accompanying publication: ["A computational tool for the analysis of 3D bending-active
structures based on the dynamic relaxation method"](https://diglib.eg.org/handle/10.2312/stag20221250).

[comment]: <> (CMake Linux:libxrandr-dev, libxinerama-dev, libxcursor-dev, libxi-dev and libeigen3-dev)
Requires: C++20,git,..
TODO:list dependencies

Currently the code has been built using: Clang 10, AppleClang 13 and MSVC 19.

![teaser_results](https://user-images.githubusercontent.com/17647952/200535216-e3cb3cba-a5c8-4ac4-bc71-5881746cc57e.png)

---
If this tool contributes to an academic publication, cite it as:
```bib
@inproceedings {manolas2022drm,
booktitle = {Smart Tools and Applications in Graphics - Eurographics Italian Chapter Conference},
editor = {Cabiddu, Daniela and Schneider, Teseo and Allegra, Dario and Catalano, Chiara Eva and Cherchi, Gianmarco and Scateni, Riccardo},
title = {{A Computational Tool for the Analysis of 3D Bending-active Structures Based on the Dynamic Relaxation Method}},
author = {Manolas, Iason and Laccone, Francesco and Cherchi, Gianmarco and Malomo, Luigi and Cignoni, Paolo},
year = {2022},
publisher = {The Eurographics Association},
ISSN = {2617-4855},
ISBN = {978-3-03868-191-5},
DOI = {10.2312/stag.20221250}
}
```
