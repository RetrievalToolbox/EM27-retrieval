#!/bin/bash

echo "Downloading required data for example.."
echo "======================================="

mkdir -p ./aux_data
cd ./aux_data

if [ -f "solar_merged_20240731_600_33300_000.out" ]; then
    echo "Solar model file exists. Skip download."
else
    echo "Downloading solar model data from JPL."
    wget https://mark4sun.jpl.nasa.gov/toon/solar/solar_merged_20240731_600_33300_000.out.gz
    gunzip solar_merged_20240731_600_33300_000.out.gz
fi

cd ..

# Install zenodo-get into venv
python3 -m venv .venv
source .venv/bin/activate
pip install zenodo-get

# Download spectroscopy from Zenodo
mkdir -p spectroscopy

zenodo_get -g '*ABSCO*.nc' -o spectroscopy -d 10.5281/zenodo.18246643

# Download the example data
zenodo_get -g 'example_data.zip' -o ./ -d 10.5281/zenodo.18246643

unzip example_data.zip # Maybe this only works on mac?