# Download required data for example
# (disclaimer: this script was written with ChatGSFC)
Write-Host "Downloading required data for example.."
Write-Host "======================================="

# Create aux_data directory if it doesn't exist
if (-not (Test-Path -Path ".\aux_data")) {
    New-Item -ItemType Directory -Path ".\aux_data" | Out-Null
}
Push-Location ".\aux_data"

# Download solar model data
if (Test-Path -Path "solar_merged_20240731_600_33300_000.out") {
    Write-Host "Solar model file exists. Skip download."
} else {
    Write-Host "Downloading solar model data from JPL."

    #$ProgressPreference = 'SilentlyContinue'
    Invoke-WebRequest -Uri "https://mark4sun.jpl.nasa.gov/toon/solar/solar_merged_20240731_600_33300_000.out.gz" `
        -OutFile "solar_merged_20240731_600_33300_000.out.gz"

    # Extract gz file
    $gzipFile = "solar_merged_20240731_600_33300_000.out.gz"
    $outputFile = "solar_merged_20240731_600_33300_000.out"
    $gzipStream = [System.IO.File]::OpenRead($gzipFile)
    $outputStream = [System.IO.File]::Create($outputFile)
    $gzipDecompressor = New-Object System.IO.Compression.GZipStream($gzipStream, [System.IO.Compression.CompressionMode]::Decompress)
    $gzipDecompressor.CopyTo($outputStream)
    $gzipDecompressor.Close()
    $outputStream.Close()
    $gzipStream.Close()
    Remove-Item $gzipFile
}

Pop-Location

# ZENODO downloads
$zenodo_record = "18246643"
$zenodo_root = "https://zenodo.org/records/${zenodo_record}/files"

# Define spectroscopy files and checksums
$spec_fnames = @(
    "CH4_04175-04350_v0.0_ABSCO.nc",
    "CH4_05850-06200_v0.0_ABSCO.nc",
    "CH4_06100-06450_v0.0_ABSCO.nc",
    "CO_04175-04350_v0.0_ABSCO.nc",
    "CO2_04780-04940_v0.0_ABSCO.nc",
    "CO2_05850-06200_v0.0_ABSCO.nc",
    "CO2_06100-06450_v0.0_ABSCO.nc",
    "H2O_04175-04350_v0.0_ABSCO.nc",
    "H2O_04780-04940_v0.0_ABSCO.nc",
    "H2O_05850-06200_v0.0_ABSCO.nc",
    "H2O_06100-06450_v0.0_ABSCO.nc",
    "H2O_07600-08200_v0.0_ABSCO.nc",
    "H2O_08250-08550_v0.0_ABSCO.nc",
    "O2_07600-08200_v0.0_ABSCO.nc"
)

$spec_md5 = @(
    "f9221c8e057009b4c85fcc350585a22a",
    "329613077f259cdff622bc8e7462a8b7",
    "5d36a39a7fa4f1afa8e0c3aa93b4a0db",
    "e501257173145a4a06f3181b7eccb2d8",
    "58f920cf86928edf2827f4ea0c66ffb0",
    "6cd6c55b6caaf1a69b111bb95f31ae13",
    "090136e55b67e43c84302c90ddfb8b87",
    "ed0ea6657c0d2cf1ddfbb23aa1300394",
    "4935e9b00e58506175c86af550727803",
    "d22216051338f4d7d434101284321f03",
    "11ed22a5f90f8d20479844427729747e",
    "ec56cf8ec0b59deb45af5ecc9305d0f4",
    "9b8a99159a864dba45be5a203265c230",
    "2407dd792aa4dc02cf0851783a11db7f"
)

# Create spectroscopy subdirectory if needed
if (-not (Test-Path -Path ".\spectroscopy")) {
    Write-Host "Creating spectroscopy subdirectory.."
    New-Item -ItemType Directory -Path ".\spectroscopy" | Out-Null
} else {
    Write-Host "spectroscopy subdirectory already exists.."
}

# Function to calculate MD5 hash
function Get-MD5Hash {
    param([string]$FilePath)
    $md5 = New-Object -TypeName System.Security.Cryptography.MD5CryptoServiceProvider
    $hash = [System.BitConverter]::ToString($md5.ComputeHash([System.IO.File]::ReadAllBytes($FilePath)))
    return $hash.Replace("-", "").ToLower()
}

# Download spectroscopy files
#$ProgressPreference = 'SilentlyContinue'
for ($i = 0; $i -lt $spec_fnames.Count; $i++) {
    $fname = $spec_fnames[$i]
    $spec_path = ".\spectroscopy\${fname}"

    if (Test-Path -Path $spec_path) {
        # Check MD5
        $computed_hash = Get-MD5Hash -FilePath $spec_path
        if ($computed_hash -eq $spec_md5[$i]) {
            Write-Host "MD5 check successful for ${fname}"
        } else {
            Write-Host "MD5 check failed!!! Manually delete and re-download ${fname}!"
        }
    } else {
        Write-Host "Downloading ${fname}"
        Invoke-WebRequest -Uri "${zenodo_root}/${fname}" -OutFile $spec_path
    }
}

# Download example data
if (Test-Path -Path ".\example_data.zip") {
    Write-Host "Example data already exists."
} else {
    Write-Host "Downloading example_data.zip"
    Invoke-WebRequest -Uri "${zenodo_root}/example_data.zip" -OutFile ".\example_data.zip"
}

# Unzip example data if needed
if (Test-Path -Path ".\example_data") {
    Write-Host "Example data directory already exists. Skipping unzipping."
} else {
    Write-Host "Unzipping example data..."
    Expand-Archive -Path ".\example_data.zip" -DestinationPath "."
    Write-Host "Done."
}

Write-Host "All done!"