#!/bin/bash
# Bash script to remap atmospheric forcing data with the mean of CDO.
# In advance the griddes file needs to be processed, e.g. by using
# cdo griddes domainfile.nc > domainfile.txt.

# Load environment 
#  ml Stages/2025  Intel/2024.2.0  ParaStationMPI/5.11.0-1 CDO/2.4.4

# Input and output directories
INPUT_DIR="/p/scratch/cslts/poll1/sim/paper/icon-cordex_forc/atm-forc/org"
OUTPUT_DIR="/p/scratch/cslts/poll1/sim/paper/icon-cordex_forc/atm-forc/"
GRIDDES_FILE="${OUTPUT_DIR}/domain_ICON-11.txt"

# Make output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"

# Loop through all files matching YYYY-MM.nc
for file in "$INPUT_DIR"/*.nc; do
#for file in "$INPUT_DIR"/195[1-9]*.nc; do
    # Extract the filename only (without path)
    fname=$(basename "$file")

    # NN interpolation using CDO
    cdo -remapnn,${GRIDDES_FILE} "$file" "$OUTPUT_DIR/$fname"

    echo "Processed: $fname"
done

