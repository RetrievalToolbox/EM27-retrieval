# EM27/SUN retrieval demo application

## Quick set-up

After cloning this repository, perform the following steps.

1. Install Julia

Download the Julia version multiplexer `JuliaUp` via

    curl -fsSL https://install.julialang.org

2. Instantiate the environment

Inside the cloned directory type following commands to install JupyterLab, needed to run the demo notebook.

    julia -e 'using Pkg(); Pkg.install("IJulia")'

Then install all required packages with (still inside the cloned directory)

    julia --project="./" -e 'using Pkg; Pkg.instantiate()'

As a quick check, running the following should produce the Earth radius without error messages:

    julia --project="./" -e 'using RetrievalToolbox; println(EARTH_RADIUS)'

3. Download auxiliary data (spectroscopy, example measurements, solar spectrum)

The demo requires additional data to be downloaded from a Zenodo repository. Run the `download_data.sh` script to automatically download the needed files:

    chmod u+x download_data.sh
    ./download_data.sh


## Launch the demo notebook

    julia # Start Julia REPL
    julia> using IJulia # Import IJulia module
    julia> IJulia.jupyterlab(dir=pwd()) # Launch JupyterLab
