"""
Convert accumulated radiation to radiative flux.

This script converts accumulated radiation values (J/m²) to instantaneous
radiative flux values (W/m²) by calculating the rate of change over the
time intervals between measurements.

Mathematical Background:
Radiative flux is the rate of energy transfer per unit area:

    F = ΔE / (A × Δt) = ΔE / Δt  [W/m²]

where:
- F: Radiative flux (W/m²) or (J/(s·m²))
- ΔE: Change in accumulated energy (J/m²)
- Δt: Time interval (s)
- A: Area (m²), which cancels out since we work with per-unit-area values

For accumulated radiation data, the flux is computed as:

    F(t) = (E_acc(t) - E_acc(t-1)) / Δt

where E_acc is the accumulated radiation at each timestep.

Notes on accumulated radiation in reanalysis products:
- ERA5 and similar reanalysis datasets provide radiation as accumulated values
- Accumulations are typically reset at the start of each forecast
- The time interval Δt should be computed from the actual time coordinates
- First timestep flux cannot be computed (requires previous accumulation)

References:
1. ECMWF IFS Documentation - Part IV: Physical Processes
2. ERA5 documentation: https://confluence.ecmwf.int/display/CKB/ERA5+documentation
"""

import argparse
import numpy as np
import netCDF4
from datetime import datetime, timedelta


def accumulated_to_flux(accumulated, time_seconds):
    """
    Convert accumulated radiation to radiative flux.

    Parameters
    ----------
    accumulated : array
        Accumulated radiation values [J/m²]
        Shape: (time, ...)
    time_seconds : array
        Time coordinates in seconds since reference
        Shape: (time,)

    Returns
    -------
    flux : array
        Radiative flux [W/m²]
        Shape: (time, ...)
        First timestep will be NaN (cannot compute without previous value)

    Notes
    -----
    The flux at timestep t is computed as:
        flux[t] = (accumulated[t] - accumulated[t-1]) / (time[t] - time[t-1])

    The first timestep will contain NaN values since there is no previous
    accumulation to compute the difference from.

    Examples
    --------
    >>> # Example with hourly data (3600 seconds interval)
    >>> accumulated = np.array([3600000, 7200000, 10800000])  # J/m²
    >>> time = np.array([0, 3600, 7200])  # seconds
    >>> flux = accumulated_to_flux(accumulated, time)
    >>> print(flux)
    [nan 1000. 1000.]  # W/m²
    """
    # Initialize flux array with same shape as accumulated
    flux = np.zeros_like(accumulated, dtype=np.float32)
    flux[0, ...] = np.nan  # First timestep cannot be computed

    # Get the number of timesteps
    n_time = accumulated.shape[0]

    # Compute flux for each timestep
    for t in range(1, n_time):
        # Calculate time difference in seconds
        dt = time_seconds[t] - time_seconds[t - 1]

        if dt <= 0:
            raise ValueError(f"Non-positive time difference at timestep {t}: dt={dt}")

        # Calculate energy difference
        de = accumulated[t, ...] - accumulated[t - 1, ...]

        # Handle potential negative values (can occur at accumulation resets)
        # If accumulated value decreases, it indicates a reset in accumulation
        # In this case, use the current accumulated value as the energy difference
        de = np.where(de < 0, accumulated[t, ...], de)

        # Compute flux: W/m² = J/m² / s
        flux[t, ...] = de / dt

    return flux


def add_radiation_flux_to_netcdf(
    filename,
    thermal_var="strd",
    solar_var="ssrd",
    time_var="time",
    thermal_flux_name="flds",
    solar_flux_name="fsds"
):
    """
    Read accumulated radiation from a netCDF file,
    calculate radiative fluxes, and write them back to the file.

    Parameters:
    -----------
    filename : str
        Path to the netCDF file
    thermal_var : str, optional
        Name of the thermal radiation variable (default: 'strd')
        Surface thermal radiation downwards (accumulated)
    solar_var : str, optional
        Name of the solar radiation variable (default: 'ssrd')
        Surface solar radiation downwards (accumulated)
    time_var : str, optional
        Name of the time variable (default: 'time')
    thermal_flux_name : str, optional
        Name for the output thermal flux variable (default: 'flds')
        Downward longwave radiation at surface
    solar_flux_name : str, optional
        Name for the output solar flux variable (default: 'fsds')
        Downward shortwave radiation at surface

    Returns:
    --------
    None
        Modifies the netCDF file in place by adding flux variables

    Raises:
    -------
    ValueError
        If flux variables already exist in the file
    KeyError
        If required variables are not found
    """

    # Open netCDF file in append mode
    print(f"Opening {filename}...")
    nc = netCDF4.Dataset(filename, "a")

    try:
        # Check if flux variables already exist
        if thermal_flux_name in nc.variables:
            nc.close()
            raise ValueError(
                f"Variable '{thermal_flux_name}' already exists in {filename}. "
                "No changes made. Delete the variable first if you want to recalculate."
            )
        if solar_flux_name in nc.variables:
            nc.close()
            raise ValueError(
                f"Variable '{solar_flux_name}' already exists in {filename}. "
                "No changes made. Delete the variable first if you want to recalculate."
            )

        # Read the required variables
        print(f"Reading {thermal_var}, {solar_var}, and {time_var}...")
        strd = nc.variables[thermal_var][:]  # Accumulated thermal radiation [J/m²]
        ssrd = nc.variables[solar_var][:]    # Accumulated solar radiation [J/m²]
        time_var_obj = nc.variables[time_var]
        time_values = time_var_obj[:]

        print(f"Data shapes - {thermal_var}: {strd.shape}, {solar_var}: {ssrd.shape}")
        print(f"Time dimension size: {len(time_values)}")

        # Convert time to seconds since first timestep
        # Handle different time units
        time_units = time_var_obj.units
        time_calendar = getattr(time_var_obj, 'calendar', 'standard')

        print(f"Time units: {time_units}")
        print(f"Time calendar: {time_calendar}")

        # Use netCDF4's num2date to convert time values to datetime objects
        time_dates = netCDF4.num2date(time_values, units=time_units, calendar=time_calendar)

        # Convert to seconds since first timestep
        time_seconds = np.array([
            (date - time_dates[0]).total_seconds() for date in time_dates
        ])

        print(f"Time range: {time_dates[0]} to {time_dates[-1]}")
        print(f"Time step (first interval): {time_seconds[1] - time_seconds[0]:.0f} seconds")

        # Calculate fluxes
        print("Calculating thermal radiation flux (longwave downward)...")
        flds = accumulated_to_flux(strd, time_seconds)

        print("Calculating solar radiation flux (shortwave downward)...")
        fsds = accumulated_to_flux(ssrd, time_seconds)

        # Create the thermal flux variable
        print(f"Creating new variable '{thermal_flux_name}'...")
        dim_names = nc.variables[thermal_var].dimensions
        flds_var = nc.createVariable(thermal_flux_name, "f4", dim_names, zlib=True, complevel=4)

        # Add attributes for thermal flux
        flds_var.units = "W m-2"
        flds_var.long_name = "Downward longwave radiation at surface"
        flds_var.standard_name = "surface_downwelling_longwave_flux_in_air"
        flds_var.description = (
            f"Calculated from accumulated {thermal_var} by differencing consecutive "
            "timesteps and dividing by the time interval"
        )
        flds_var.missing_value = np.nan
        flds_var[:] = flds

        print(f"Variable '{thermal_flux_name}' created and written successfully!")

        # Create the solar flux variable
        print(f"Creating new variable '{solar_flux_name}'...")
        dim_names = nc.variables[solar_var].dimensions
        fsds_var = nc.createVariable(solar_flux_name, "f4", dim_names, zlib=True, complevel=4)

        # Add attributes for solar flux
        fsds_var.units = "W m-2"
        fsds_var.long_name = "Downward shortwave radiation at surface"
        fsds_var.standard_name = "surface_downwelling_shortwave_flux_in_air"
        fsds_var.description = (
            f"Calculated from accumulated {solar_var} by differencing consecutive "
            "timesteps and dividing by the time interval"
        )
        fsds_var.missing_value = np.nan
        fsds_var[:] = fsds

        print(f"Variable '{solar_flux_name}' created and written successfully!")

        # Print some statistics (excluding first timestep with NaN)
        print("\nStatistics:")
        print(f"\nThermal radiation flux ({thermal_flux_name}):")
        print(f"  Range: {np.nanmin(flds):.2f} to {np.nanmax(flds):.2f} W/m²")
        print(f"  Mean: {np.nanmean(flds):.2f} W/m²")
        print(f"  Std: {np.nanstd(flds):.2f} W/m²")

        print(f"\nSolar radiation flux ({solar_flux_name}):")
        print(f"  Range: {np.nanmin(fsds):.2f} to {np.nanmax(fsds):.2f} W/m²")
        print(f"  Mean: {np.nanmean(fsds):.2f} W/m²")
        print(f"  Std: {np.nanstd(fsds):.2f} W/m²")

        # Check for negative values (excluding NaN)
        n_negative_thermal = np.sum(flds[~np.isnan(flds)] < 0)
        n_negative_solar = np.sum(fsds[~np.isnan(fsds)] < 0)

        if n_negative_thermal > 0:
            print(f"\nWarning: {n_negative_thermal} negative values in thermal flux")
        if n_negative_solar > 0:
            print(f"\nWarning: {n_negative_solar} negative values in solar flux")

        print(f"\nNote: First timestep contains NaN values (no previous accumulation)")

    except KeyError as e:
        print(f"Error: Required variable not found in netCDF file: {e}")
        print(f"Available variables: {list(nc.variables.keys())}")
        nc.close()
        raise

    except ValueError as e:
        print(f"Error: {e}")
        raise

    except Exception as e:
        print(f"Unexpected error occurred: {e}")
        nc.close()
        raise

    else:
        # Only executes if no exception was raised
        nc.close()
        print(f"\nFile {filename} closed successfully.")


if __name__ == "__main__":
    # Set up argument parser
    parser = argparse.ArgumentParser(
        description="Add radiative flux variables to a netCDF file based on accumulated radiation."
    )
    parser.add_argument(
        "filename",
        type=str,
        help="Path to the netCDF file containing accumulated radiation variables"
    )
    parser.add_argument(
        "--thermal-var",
        type=str,
        default="strd",
        help="Name of the accumulated thermal radiation variable (default: strd)"
    )
    parser.add_argument(
        "--solar-var",
        type=str,
        default="ssrd",
        help="Name of the accumulated solar radiation variable (default: ssrd)"
    )
    parser.add_argument(
        "--time-var",
        type=str,
        default="time",
        help="Name of the time variable (default: time)"
    )
    parser.add_argument(
        "--thermal-flux-name",
        type=str,
        default="flds",
        help="Name for the output thermal flux variable (default: flds)"
    )
    parser.add_argument(
        "--solar-flux-name",
        type=str,
        default="fsds",
        help="Name for the output solar flux variable (default: fsds)"
    )

    # Parse command-line arguments
    args = parser.parse_args()

    # Process the file
    add_radiation_flux_to_netcdf(
        args.filename,
        thermal_var=args.thermal_var,
        solar_var=args.solar_var,
        time_var=args.time_var,
        thermal_flux_name=args.thermal_flux_name,
        solar_flux_name=args.solar_flux_name
    )
