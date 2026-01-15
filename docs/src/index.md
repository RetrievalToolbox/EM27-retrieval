# EM27/SUN Retrieval Algorithm at NASA Goddard Space Flight Center


## Introduction

This software implements a flexible retrieval algorithm for use with spectra acquired by EM27/SUN measurement devices. The provided package includes a batch script that can be used to process a series of measurements, and a demonstration Jupyter notebook that might be helpful to understand the steps of the processing which are more hidden in the batch script.

## Requirements

This package uses `G3RT.jl` (Goddard Greenhouse Gas Retrieval Toolkit) to implement the retrieval algorithm, which must be pre-installed.

!!! warning "TODO"
    Add working instructions as to how to obtain G3RT.jl from the NASA GitHub repo.

Currently, the algorithm is set up to process spectra with have been pre-processed with the `PROFFAST`[^PROFFAST] tools, and prior atmosphere profiles in the `*.map` format, as provided by `Ginput`[^Ginput], the GGG meteorology and a priori VMR preprocessor.

Embedding this software into a processing workflow for EM27/SUN measurements would entail roughly the following steps:

1. Acquire measurements with the EM27/SUN instrument, and optionally other meteorological data such as surface pressure
2. Run the PROFFAST pre-processor on a set of interferograms that are stored in a directory
3. Obtain the TCCON-Ginput `*.map` files through the TCCON service or by running the software yourself (which needs the appropriate meteorological model data)
4. Execute the batch script

[^PROFFAST]: [PROFFAST at Karlsruhe Institute for Technology](https://www.imk-asf.kit.edu/english/3225.php)
[^Ginput]: [https://github.com/TCCON/py-ginput](https://github.com/TCCON/py-ginput)