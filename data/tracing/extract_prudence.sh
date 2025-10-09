#!/bin/bash

# This script merges and extracts the PRUDENCE regions from tracing data as given in this publication:
# Mayer, Amelie, and Volkmar Wirth. "Lagrangian characterization of heat waves: The perspective matters." EGUsphere 2024 (2024): 1-29.
# The data is available at: https://doi.org/10.5281/zenodo.14679142
# This script uses the files contained in "absolutes_YYYY.tar"

# PRUDENCE region definitions (lon1,lon2,lat1,lat2)
declare -A regions=(
    ["BI"]="-10,2,50,59"      # British Isles
    ["IP"]="-10,3,36,44"      # Iberian Peninsula
    ["FR"]="-5,5,44,50"       # France
    ["ME"]="2,16,48,55"       # Mid-Europe
    ["SC"]="5,30,55,70"       # Scandinavia
    ["AL"]="5,15,44,48"       # Alps
    ["MD"]="3,25,36,44"       # Mediterranean
    ["EA"]="16,30,44,55"      # Eastern Europe
)

export -A regions

process_year() {
    year=$1
    echo "Processing year $year"
  
    # Merge summer months (JJA)
    cdo mergetime ${year}/temperature_diagnostics_${year}-0[6-8]-??.nc temp_${year}_JJA.nc
    
    # Extract each PRUDENCE region
    for region in "${!regions[@]}"; do
        bbox=${regions[$region]}
        cdo sellonlatbox,$bbox temp_${year}_JJA.nc temperature_diagnostics_${region}_${year}.nc
    done
    
    # Clean up temporary file
    rm temp_${year}_JJA.nc
}


for year in {2010..2022}; do
    process_year $year &
done
wait  # Wait for all background jobs to complete
