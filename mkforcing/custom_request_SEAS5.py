## Download for SEAS5 forecast variables
# https://cds.climate.copernicus.eu/datasets/seasonal-original-single-levels?tab=download

# Forecast for 7 months

# import cdsapi

dataset = "seasonal-original-single-levels"
request = {
    "originating_centre": "ecmwf",
    "system": "51",
    "variable": [
        "mean_sea_level_pressure",  # convert to surface pressure, use elevation (hypsometric formula), surface geopotential height (orography)
        "orography",  # used to convert mslp to sp
        "surface_thermal_radiation_downwards",  # Unit conversion from accumulated value [J/m2] to mean rate [W/m2]
        "surface_solar_radiation_downwards",  # Unit conversion from accumulated value [J/m2] to mean rate [W/m2]
        "total_precipitation",
        "10m_u_component_of_wind",
        "10m_v_component_of_wind",
        "2m_temperature",
        "2m_dewpoint_temperature",
    ],
    "year": ["2025"],
    "month": ["09"],
    "day": ["01"],
    "leadtime_hour": [str(h) for h in range(6, 5161, 6)],
    "data_format": "netcdf",
    "area": [51, 6, 50, 7]      # Selhausen
    # "area": [74, -42, 20, 69] # Europe
}

# client = cdsapi.Client()
# client.retrieve(dataset, request).download()
