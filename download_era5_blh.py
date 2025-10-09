# Download script for ERA5 BLH data via the CDS API.
# Further, we postprocessed this data by extracting daily maxima
# and merging into a single file. 

import cdsapi

client = cdsapi.Client()

for year in [2000, 2024]:
    dataset = "reanalysis-era5-single-levels"
    request = {
        "product_type": ["reanalysis"],
        "variable": ["boundary_layer_height"],
        "year": [str(year)],
        "month": ["06", "07", "08"],
        "day": [
            "01", "02", "03",
            "04", "05", "06",
            "07", "08", "09",
            "10", "11", "12",
            "13", "14", "15",
            "16", "17", "18",
            "19", "20", "21",
            "22", "23", "24",
            "25", "26", "27",
            "28", "29", "30",
            "31"
        ],
        "time": [
            "08:00", "09:00", "10:00",
            "11:00", "12:00", "13:00",
            "14:00", "15:00", "16:00",
            "17:00", "18:00"
        ],
        "data_format": "netcdf",
        "download_format": "unarchived",
        "area": [75, -25, 30, 60]
    }
    target = f"blh_era5_{year}.nc"
    client.retrieve(dataset, request, target)#.download()
