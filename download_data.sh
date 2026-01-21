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

## ZENODO downloads
zenodo_record=18246643
zenodo_root=https://zenodo.org/records/${zenodo_record}/files

# Download spectroscopy from Zenodo
declare -a spec_fnames
spec_fnames[0]="CH4_04175-04350_v0.0_ABSCO.nc"
spec_fnames[1]="CH4_05850-06200_v0.0_ABSCO.nc"
spec_fnames[2]="CH4_06100-06450_v0.0_ABSCO.nc"
spec_fnames[3]="CO_04175-04350_v0.0_ABSCO.nc"
spec_fnames[4]="CO2_04780-04940_v0.0_ABSCO.nc"
spec_fnames[5]="CO2_05850-06200_v0.0_ABSCO.nc"
spec_fnames[6]="CO2_06100-06450_v0.0_ABSCO.nc"
spec_fnames[7]="H2O_04175-04350_v0.0_ABSCO.nc"
spec_fnames[8]="H2O_04780-04940_v0.0_ABSCO.nc"
spec_fnames[9]="H2O_05850-06200_v0.0_ABSCO.nc"
spec_fnames[10]="H2O_06100-06450_v0.0_ABSCO.nc"
spec_fnames[11]="H2O_07600-08200_v0.0_ABSCO.nc"
spec_fnames[12]="H2O_08250-08550_v0.0_ABSCO.nc"
spec_fnames[13]="O2_07600-08200_v0.0_ABSCO.nc"

declare -a spec_md5
spec_md5[0]="f9221c8e057009b4c85fcc350585a22a"
spec_md5[1]="329613077f259cdff622bc8e7462a8b7"
spec_md5[2]="5d36a39a7fa4f1afa8e0c3aa93b4a0db"
spec_md5[3]="e501257173145a4a06f3181b7eccb2d8"
spec_md5[4]="58f920cf86928edf2827f4ea0c66ffb0"
spec_md5[5]="6cd6c55b6caaf1a69b111bb95f31ae13"
spec_md5[6]="090136e55b67e43c84302c90ddfb8b87"
spec_md5[7]="ed0ea6657c0d2cf1ddfbb23aa1300394"
spec_md5[8]="4935e9b00e58506175c86af550727803"
spec_md5[9]="d22216051338f4d7d434101284321f03"
spec_md5[10]="11ed22a5f90f8d20479844427729747e"
spec_md5[11]="ec56cf8ec0b59deb45af5ecc9305d0f4"
spec_md5[12]="9b8a99159a864dba45be5a203265c230"
spec_md5[13]="2407dd792aa4dc02cf0851783a11db7f"


fcount=0
for fname in "${spec_fnames[@]}"; do

    # Check if file is there
    if [ -f "spectroscopy/${fname}" ]; then

        # Check MD5
        if echo "$spec_md[$fcount]" spectroscopy/${fname} | md5sum --status -c -; then
            echo "MD5 check successful for ${fname}"
        else
            echo "MD5 check failed!!! Manually delete and re-download ${fname}!"
        fi

    else
        echo "Downloading ${fname}"
        curl --remove-on-error -o spectroscopy/${fname} ${zenodo_root}/${fname}
    fi

    fcount=$((fcount + 1))
done


# Download the example data
if [ -f "example_data.zip" ]; then
    echo "Example data already exists."
else
    curl ${zenodo_root}/example_data.zip -o example_data.zip
fi

# Unzip example data if unzipped directory does not exist
if [ -d "example_data" ]; then
    echo "Example data directory already exists. Skipping unzipping."
else
    echo "Unzipping example data..."
    unzip -q example_data.zip
    echo "Done."
fi

echo "All done!"