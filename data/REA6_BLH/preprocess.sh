#!/bin/bash

# Process data from disk:
for year in {2014..2018}; do
  echo "Combine and compute daily maxima for year $year"
  #Merge and rechunk netCDF
  cdo -b f32 -z zip2 -k auto -w mergetime "${year}/blh_0???_bri.nc" "blh_${year}_bri.nc"
  cdo -b f32 -z zip2 -k auto -w mergetime "${year}/blh_0???_bri_int_alt.nc" "blh_${year}_bri_int_alt.nc"

  #Extract daily maxima
  cdo -b f32 -z zip2 -k auto daymax "blh_${year}_bri.nc" "blh_dmax_${year}_bri.nc"
  cdo -b f32 -z zip2 -k auto daymax "blh_${year}_bri_int_alt.nc" "blh_dmax_${year}_bri_int_alt.nc"
done