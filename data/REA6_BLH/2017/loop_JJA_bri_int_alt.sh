#!/bin/bash
year=${1}
echo "Starting estimation procedure for ${year} @ $(date +"%T") using the Bulk Richardson Method"
echo


# Array of months
month_names=("June" "July" "August")

# Loop through each month
for month_name in "${month_names[@]}"; do
    echo "Month: $month_name"


    # Determine the number of days in the month
    case $month_name in
        "June")
            days_in_month=30
            month=06
            ;;
        "July")
            days_in_month=31
            month=07
            ;;
        "August")
            days_in_month=31
            month=08
            ;;
    esac
    

    # Loop through all remaining days in the month
    for ((day=1; day<=$days_in_month; day++)); do
        echo "  Day: $day"
        printf -v out '%02d' "$day"


        # Retrieve data from archive
        ./extract_tars.sh $year $month $out
        if [ $? -eq 0 ]; then
          echo "    - Extraction finished @ $(date +"%T")"
        else
          echo "    - Extraction failed"
          continue
        fi


        # Fit the boundary layer height
        TIMEOUT=600
        MAX_RETRIES=2
        RETRY_COUNT=0
        while [ $RETRY_COUNT -lt $MAX_RETRIES ]; do
          timeout $TIMEOUT python3 ../bri_fit_int_alt.py $year $month $out

          if [ $? -eq 0 ]; then
            echo "    - Estimation finished @ $(date +"%T")"
            break
          else
            echo "    - Estimation failed (${RETRY_COUNT+1} / ${MAX_RETRIES})"
            ((RETRY_COUNT++))
          fi
        done
      

        # Do clean-up
        if [ -r "./blhs_merged.nc" ]
        then
          rm "blhs.nc"
          mv "blhs_merged.nc" "blhs.nc"
        fi
        
        rm ./PP/PP.3D.${year}*.nc
        rm ./T/T.3D.${year}*.nc
        rm ./QV/QV.3D.${year}*.nc
        rm ./U/U.3D.${year}*.nc
        rm ./V/V.3D.${year}*.nc
        echo "    - Clean-up done"
    done

    mv "blhs.nc" "blh_${month}${year: -2}_bri_int_alt.nc"

    echo # Add a newline between months
done
