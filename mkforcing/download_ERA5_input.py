#!/usr/bin/env python3
"""
Download ERA5 reanalysis data from Copernicus Climate Data Store (CDS).

This script downloads ERA5 single-level reanalysis data for a specified month
and year. The default data includes surface pressure, radiation fluxes, and
precipitation at hourly resolution.

Requirements:
    - cdsapi library (pip install cdsapi)
    - CDS API credentials configured in ~/.cdsapirc

Usage:
    python download_ERA5_input.py --year <year> --month <month> --dirout <output_directory>
    python download_ERA5_input.py --year 2017 --month 7 --dirout ./output
    python download_ERA5_input.py --year 2017 --month 7 --dirout ./output --request custom_request.py
    python download_ERA5_input.py --help

Note:
    CDS API credentials must be configured before use.
    See: https://cds.climate.copernicus.eu/api-how-to
"""
import argparse
import calendar
import cdsapi
import sys
import os

def generate_days(year, month):
    """Get the number of days in a given month and year.

    Args:
        year (int): Year
        month (int): Month (1-12)

    Returns:
        list: List of day numbers for the month
    """
    # Get the number of days in the given month
    num_days = calendar.monthrange(year, month)[1]

    # Generate the list of days as integers
    days = [day for day in range(1, num_days + 1)]

    return days

def generate_datarequest(year, monthstr, days,
                        dataset="reanalysis-era5-single-levels",
                        request=None,
                        target=None):
    """Generate and execute ERA5 data download request.

    "ERA5 hourly data on single levels from 1940 to present":
    https://cds.climate.copernicus.eu/datasets/reanalysis-era5-single-levels?tab=overview

    Args:
        year (int): Year to download
        monthstr (str): Month as zero-padded string (e.g., '07')
        days (list): List of days in the month
        dataset (str, optional): CDS dataset name. Defaults to 'reanalysis-era5-single-levels'.
        request (dict, optional): Custom CDS request dictionary. If None, uses default request.
        target (str, optional): Output filename. If None, uses 'download_era5_YYYY_MM.zip'.

    Returns:
        str: Path to downloaded file
    """

    # active download client for climate data service (cds)
    client = cdsapi.Client()

    # Default request if not provided
    if request is None:
        request = {
            "product_type": ["reanalysis"],
            "variable": [
                "surface_pressure",
                "mean_surface_downward_long_wave_radiation_flux",
                "mean_surface_downward_short_wave_radiation_flux",
                "mean_total_precipitation_rate"
            ],
            "year": [str(year)],
            "month": [monthstr],
            "day": days,
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
            "area": [74, -42, 20, 69]
        }
    else:
        # Adapt year, month and day to input values
        request["year"] = [str(year)]
        request["month"] = [monthstr]
        request["day"] = days

    # Default filename if not provided
    if target is None:
        target = f'download_era5_{year}_{monthstr}.zip'

    # Get the data from cds
    client.retrieve(dataset, request, target)

    return target


if __name__ == "__main__":

    # Set up argument parser
    parser = argparse.ArgumentParser(
        description="Download ERA5 reanalysis data from Copernicus Climate Data Store (CDS)."
    )
    parser.add_argument(
        "--year",
        type=int,
        required=True,
        help="Year to download (e.g., 2017)"
    )
    parser.add_argument(
        "--month",
        type=int,
        required=True,
        help="Month to download (1-12)"
    )
    parser.add_argument(
        "--dirout",
        type=str,
        required=True,
        help="Output directory path"
    )
    parser.add_argument(
        "--request",
        type=str,
        required=False,
        help="Path to Python file defining custom 'request' variable"
    )
    parser.add_argument(
        "--dataset",
        type=str,
        required=False,
        default="reanalysis-era5-single-levels",
        help="CDS dataset name (default: reanalysis-era5-single-levels)"
    )

    # Parse command-line arguments
    args = parser.parse_args()
    year = args.year
    month = args.month
    dirout = args.dirout

    # Load custom request if provided
    custom_request = None
    custom_dataset = args.dataset
    if args.request:
        import importlib.util
        spec = importlib.util.spec_from_file_location("custom_request_module", args.request)
        custom_module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(custom_module)
        if hasattr(custom_module, 'request'):
            custom_request = custom_module.request
            print(f"Loaded custom request from: {args.request}")
        else:
            print(f"Warning: No 'request' variable found in {args.request}, using default")
        if hasattr(custom_module, 'dataset'):
            custom_dataset = custom_module.dataset
            print(f"Loaded custom dataset from: {args.request}")

    # Ensure the output directory exists, if not, create it
    if not os.path.exists(dirout):
        os.makedirs(dirout)

    # change to output directory
    os.chdir(dirout)

    # Format the month with a leading zero if needed
    monthstr = f"{month:02d}"

    # Get the list of days for the request
    days = generate_days(year, month)

    print(f"Downloading ERA5 data for {year}-{monthstr}")
    print(f"Dataset: {custom_dataset}")
    print(f"Output directory: {os.getcwd()}")

    # Execute download request
    target = generate_datarequest(year, monthstr, days, dataset=custom_dataset, request=custom_request)
    print(f"Download complete: {os.path.abspath(target)}")
