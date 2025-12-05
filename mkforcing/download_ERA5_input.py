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
    python download_ERA5_input.py <year> <month> <output_directory>
    python download_ERA5_input.py 2017 7 ./output
    python download_ERA5_input.py --help

Note:
    CDS API credentials must be configured before use.
    See: https://cds.climate.copernicus.eu/api-how-to
"""
import calendar
import cdsapi
import sys
import os
import tempfile


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


def detect_file_type(filepath):
    """Detect if downloaded file is NetCDF, GRIB or ZIP format.

    Args:
        filepath (str): Path to the downloaded file

    Returns:
        str: File extension ('.nc' for NetCDF, '.zip' for ZIP, '.grib' for GRIB)

    Raises:
        ValueError: If file format is not recognized as NetCDF, GRIB, or ZIP
    """
    # Read file magic bytes
    with open(filepath, 'rb') as f:
        magic = f.read(8)

    # ZIP files start with 'PK' (0x504B)
    if magic[:2] == b'PK':
        return '.zip'

    # NetCDF files start with 'CDF' (0x43444601 or 0x43444602) or HDF5 signature
    if magic[:3] == b'CDF' or magic[:4] == b'\x89HDF':
        return '.nc'

    # GRIB files start with 'GRIB'
    if magic[:4] == b'GRIB':
        return '.grib'

    # If we reach here, the file format is not recognized
    magic_hex = magic.hex()
    raise ValueError(
        f"Unrecognized file format for '{filepath}'. "
        f"Magic bytes: {magic_hex}. "
        f"Expected NetCDF (CDF/HDF5), GRIB, or ZIP (PK) format."
    )


def generate_datarequest(year, monthstr, days,
                         dataset="reanalysis-era5-single-levels",
                         request=None,
                         target=None):
    """Generate and execute ERA5 data download request.

    Args:
        year (int): Year to download
        monthstr (str): Month as zero-padded string (e.g., '07')
        days (list): List of days in the month
        dataset (str, optional): CDS dataset name. Defaults to 'reanalysis-era5-single-levels'.
        request (dict, optional): Custom CDS request dictionary. If None, uses default request.
        target (str, optional): Output filename. If None, auto-detects extension based on downloaded file type.

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

    # Temporary filename w/o extension if not provided
    auto_detect_extension = target is None
    if auto_detect_extension:
        # Create a temporary file for download
        temp_fd, target = tempfile.mkstemp(
            prefix=f'download_era5_{year}_{monthstr}',
            dir='.')
        os.close(temp_fd)  # Close the file descriptor

    # Get the data from cds
    client.retrieve(dataset, request, target)

    # If target was not provided, detect the file type after download
    if auto_detect_extension:
        # Detect the actual file type
        extension = detect_file_type(target)

        # Rename to final target with correct extension
        final_target = f'{target}{extension}'
        os.rename(target, final_target)
        target = final_target

    return target


if __name__ == "__main__":

    # Get the year and month from command-line arguments
    year = int(sys.argv[1])
    month = int(sys.argv[2])
    dirout = sys.argv[3]

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
    print(f"Output directory: {os.getcwd()}")

    # Execute download request
    target = generate_datarequest(year, monthstr, days)
    print(f"Download complete: {os.path.abspath(target)}")
