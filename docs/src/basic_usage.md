# Basic Usage

## Performing the Required Pre-processing Steps

The following sections explain the steps needed to produce the required input data for the retrieval algorithm. It is assumed that users have successfully produced measurements with EM27/SUN device(s), as well as other neccessary meteorological measurements, with surface pressure being the one that is strongly needed for accurate retrievals.

### PROFFAST Preprocessor

As of now, this algorithm does not process the measured interferograms, but the spectra obtained from them. Specifically, the algorithm reads files that are produced by the PROFFAST preprocessor program, which is freely available from the Karlsruhe Institute of Technology[^PROFFAST]. We recommend to download the PROFFAST code, which broadly consists of three modules - the preprocessor, a tool to calculate scene-dependent molecular absorption coefficients, as well as the algorithm to retrieve the gas concentrations from pre-processed spectra.

#### Downloading PROFFAST

PROFFAST is available at the following link [https://www.imk-asf.kit.edu/english/3225.php](https://www.imk-asf.kit.edu/english/3225.php). We strongly suggest to download both the PROFFAST algorithm, as well as `PROFFASTpylot` available on the same site, which is a Python front-end that provides a more streamlined interface to running the various PROFFAST modules.

#### Compiling on Linux


#### Compiling on Windows


#### Compiling on Mac

### TCCON MAP Atmosphere Data

## Calling the Batch Script

The provided batch script will go through all pre-processed spectra within a directory, and perform retrievals on those spectra based on the spectral window configuration that was also supplied via command line arguments. The batch script is called as follows on a Linux or Mac system and a `BASH`-type shell:

```bash
SPEC_FOLDER=./example_spec
MAP_FOLDER=./example_map
WINDOWS_CONFIG=windows.yml
WINDOWS="2 3"
OUTPUT="example_results"

julia EM27-run.jl \
    --folder_spectra ${SPEC_FOLDER} \
    --folder_map ${SPEC_FOLDER} \
    --windows_config ${WINDOWS_CONFIG} \
    --windows ${WINDOWS} \
    --output_basename ${OUTPUT}
```

### Directories Containing Spectra and map files

The directory that is passed via the command line argument `--folder_spectra` must have the pre-processed spectra in them, and they must end in `*.BIN`. The script does not search sub-directories contained within. The same is true for the command line argument `--folder_map`, which searches the supplied directory for `*.map` files, and does not look into sub-directories. Ideally, users will not rename files that the PROFFAST pre-processor produces, which are automatically named `*??.BIN`, where the two final characters of the filename (without extension) describe the detector used during the acquisition of the interferogram (e.g., `SN` or `SM`).

### Spectral Window Configuration Selection

The window definition file is explained in more detail in [Spectral Windows File](@ref). In it, users can add or modify the spectral window limits or the gases (and their respective spectroscopy data) to be retrieved from those spectral windows. Within the spectral windows file, users can define several retrieval window configurations that are assigned by integer values. For example, in the default `windows.yml` file that is included in this software package, the window configuration `2` corresponds to the configuration to retrieve O$_2$ from the ~7800 cm$^{-1}$ wavenumber region.

The batch script accepts either a single number, or a sequence (not needed to be ordered) of numbers separated by spaces, which the script will process in the supplied order. Note that the script will exit with an error if duplicates of a window configuration number are supplied. A valid entry would be, for example:

```bash
--windows 2 3 4
```

whereas this would cause the batch script to exit with an error:

```bash
--windows 2 3 3 4
```

### Output File Basename

The argument `--output_basename` determines the filename prefix of the output file in which the result table will be written. The spectral window configuration index will be appended to the output filename to create a file for each. For example,
```bash
--output_basename my_result
--windows 3 4
```

will produce two different files: `my_result_3.csv` and `my_result_4.csv`.
!!! warning
    Existing files will be overwritten without further warning or acknowledgement!

This command line argument is optional, and the default value is `output`.

## Analysis and Plots




[^PROFFAST]: [PROFFAST at Karlsruhe Institute for Technology](https://www.imk-asf.kit.edu/english/3225.php)