# Algorithm Details

This section lists the details of the GSFC EM27/SUN retrieval algorithm, in particular the batch script, and how those can be adapted to user needs.


## Spectral Windows Definitions

### Spectral Windows File

When invoking the batch script, such as

```bash
EM27-run.jl \
    --folder_spectra ./example_data \
    --folder_map ./example_data \
    --windows 1 \
    --windows_config windows.yml
    --output_basename example
```
the commandline argument `--windows-config` is used to point to a YAML file that contains the definitions of the spectral windows that are known to the algorithm at run time. This window configuration file must have following structure, below is an example.

```yaml
1:
    name: "H2O"
    wavenumber_start: 8353.4
    wavenumber_end: 8463.1
    wavenumber_reference: 8400.0
    detector: "SN"
    spectroscopy:
        - name: "H2O"
          ABSCO: "./HITRAN_H2O_8250_8550.h5"
          unit: "1"

```

Each spectral window section must be a unique integer-valued number. Within the section the following variables must be defined: 
- `name`, which can be any type of (non-unique) identifier
- `wavenumber_start` and `wavenumber_end`, the start and end points of the spectral window in wavenumber units
- `wavenumber_reference` is a spectral reference point used in various state vector elements such as the slope of the solar scale factor
-  `detector` is either `SN` or `SM` (for now) and determines which spectrum file is being associated with this spectral window due to which detector it was recorded on (e.g. the methane band near 4250 cm$^{-1}$ is stored in a `*SM.bin` file, where as the CO$_2$ band near 6300 cm$^{-1}$ is stored in a `*SN.bin` file

The `spectroscopy` entry itself requires sub-sections, one sub-section for each gas that shall be retrieved. These elements are required for each `spectroscopy` entry:

- `name`, the name associated with the gas
- `ABSCO`, the path to the spectroscopy table which will be attached to the gas
- `unit`, a `Unitful`-type unit that will be parsed to indicate the unit of the volume mixing ratio for this gas

Note that double-quotes for strings are not strictly necessary with the `YAML.jl` package, however it makes the intention very clear that the entry should be a string-type. More importantly **ensure that the `unit` entry is always double-quoted to guarantee that it is parsed as a string**. The `unit` field can be any dimensionless `Unitful` type, and the algorithm will store and process the gas volume mixing ratio in that unit. Examples are: `"1"` to keep the volume mixing ratio simply in parts; `"ppm"` to use parts-per-million, `"ppb"` for parts-per-bilion and so forth.

Users can add new sections, or modify existing sections in this file and set up a bigger library of different spectral windows with various combinations of spectral window ranges and spectroscopy data. The batch program itself will only select the spectral windows that are passed via the `--windows` command line argument, and ignores all other spectral windows that are found in the file.

### Internal Splitting of Measurements

Internally, the algorithm reads in a full spectrum file, which can contain up to ~35,000 spectral points, depending on the detector. To make the management of measurement data easier to handle, the spectra are cut down to the size of the user-requested spectral windows, include an appropriate buffer on either side (default $\Delta\nu = 25.0\;\text{cm}^{-1}$). This is done by taking a copy of the sub-selected spectrum which is then assigned to the respective spectral window and used for the rest of the algorithm. Since those new spectrum slices are copies of the original data, users can specify several spectral windows with overlapping ranges - each of the spectral windows are processed independently. This functionality could be used to, for example, perform multiple retrievals of the same spectral window and gases, but using different spectroscopy - in a single execution of the script. Note, however, that the batch script loads the spectroscopy data for all used spectral windows at the start of the process, so every added spectral window will increase the runtime memory footprint.

As a side effect, the ISRF tables, which at the moment must be of the same size as the number of measurement points, will be smaller for each retrieval window which in turn reduces the overall memory footprint.