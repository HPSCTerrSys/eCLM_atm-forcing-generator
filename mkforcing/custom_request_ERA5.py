# Custom ERA5 request, when all information should be downloaded from ERA5
# https://cds.climate.copernicus.eu/datasets/reanalysis-era5-single-levels?tab=download

# import cdsapi

dataset = "reanalysis-era5-single-levels"
request = {
    "product_type": ["reanalysis"],
    "variable": [
        "surface_pressure",
        "mean_surface_downward_long_wave_radiation_flux",
        "mean_surface_downward_short_wave_radiation_flux",
        "mean_total_precipitation_rate",
        "10m_u_component_of_wind",
        "10m_v_component_of_wind",
        "2m_temperature",
        "2m_dewpoint_temperature",
    ],
    "time": [
        "00:00", "01:00", "02:00",
        "03:00", "04:00", "05:00",
        "06:00", "07:00", "08:00",
        "09:00", "10:00", "11:00",
        "12:00", "13:00", "14:00",
        "15:00", "16:00", "17:00",
        "18:00", "19:00", "20:00",
        "21:00", "22:00", "23:00"
    ],
    "data_format": "netcdf",
    "download_format": "unarchived",
    "area": [50.870906, 6.4421445, 50.870906, 6.4421445]      # Selhausen
    # "area": [74, -42, 20, 69] # Europe
}

# client = cdsapi.Client()
# client.retrieve(dataset, request).download()
