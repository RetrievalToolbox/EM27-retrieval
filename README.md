# EM27/SUN retrieval demo application

## Quick set-up

After cloning this repository, perform the following steps.

1. Install Julia

Download the Julia version multiplexer `JuliaUp` via

    curl -fsSL https://install.julialang.org | sh

The install will prompt users to re-source their shell, or re-connect to the machine for a fresh login. This is necessary for `JuliaUp` to work correctly for the next steps!

2. Instantiate the environment

Inside the cloned directory type following commands to install all required packages:

    julia --project="./" -e 'using Pkg; Pkg.instantiate()'

This step can take several minutes to finish. Julia will install all required packages for this demo, which includes IJulia, the Julia Jupyter kernel, including all dependencies.
As a quick check, running the following should produce the Earth radius without error messages:

    julia --project="./" -e 'using RetrievalToolbox; println(EARTH_RADIUS)'

3. Download auxiliary data (spectroscopy, example measurements, solar spectrum)

The demo requires additional data to be downloaded from a Zenodo repository. Run the `download_data.sh` script to automatically download the needed files:

    chmod u+x download_data.sh
    ./download_data.sh

This step downloads the solar model as well as the spectroscopy tables needed for the various spectral windows features in this demo. The files are downloaded from Zenodo, for a total of ~6.5 GB.

## Launch the demo notebook

    julia # Start Julia REPL
    julia> using IJulia # Import IJulia module
    julia> IJulia.jupyterlab(dir=pwd()) # Launch JupyterLab
