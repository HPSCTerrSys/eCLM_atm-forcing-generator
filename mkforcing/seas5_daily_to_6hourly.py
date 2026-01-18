#!/usr/bin/env python3
"""
Convert SEAS5 constant and daily variables to 6-hourly resolution and merge
with existing 6-hourly file.

This script:
1. Adds constant `z` (orography) to the 6-hourly file, broadcast to all time steps
2. Converts daily `tp` (total precipitation) to 6-hourly by dividing by 4
3. Converts daily `strd` and `ssrd` (radiation) to 6-hourly:
   - Zero in first/last 6-hour intervals ("night")
   - Half the daily value in the middle two intervals ("day")

Usage:
    python seas5_daily_to_6hourly.py --const <const_file> --daily <daily_file> --hourly <6h_file> [--output <output_file>]

If --output is not specified, the 6-hourly file will be modified in place.
"""

import argparse
import os
import sys

import numpy as np
import xarray as xr


def forecast_period_to_hours(forecast_period):
    """
    Convert forecast_period values to hours.

    xarray may store forecast_period as timedelta64[ns] (nanoseconds),
    but we need hours for our calculations.

    Parameters
    ----------
    forecast_period : array-like
        Forecast period values (may be in nanoseconds or hours)

    Returns
    -------
    np.ndarray
        Forecast period values in hours
    """
    values = np.asarray(forecast_period)

    # Check if values are in nanoseconds (very large numbers)
    # 1 hour = 3600 seconds = 3.6e12 nanoseconds
    if values.dtype.kind == "m":  # timedelta type
        # Convert timedelta64 to hours
        return values.astype("timedelta64[h]").astype(float)
    elif np.max(values) > 1e9:
        # Likely nanoseconds, convert to hours
        return values / (3600 * 1e9)
    else:
        # Assume already in hours
        return values.astype(float)


def expand_constant_to_6hourly(ds_const, n_timesteps):
    """
    Expand constant variable z to 6-hourly resolution.

    The constant z has forecast_period=0 (single time), we need to broadcast
    it to all 6-hourly time steps.

    Parameters
    ----------
    ds_const : xarray.Dataset
        Dataset containing constant variables (z)
    n_timesteps : int
        Number of 6-hourly time steps to expand to

    Returns
    -------
    xarray.DataArray
        z variable expanded to n_timesteps
    """
    z = ds_const["z"]

    # Remove the forecast_period dimension (size 1) and we'll broadcast later
    z_squeezed = z.squeeze("forecast_period", drop=True)

    # Expand along a new forecast_period dimension
    z_expanded = z_squeezed.expand_dims("forecast_period", axis=2)
    z_expanded = z_expanded.broadcast_to(
        number=z_squeezed.sizes.get("number", z_squeezed.shape[0]),
        forecast_reference_time=z_squeezed.sizes.get("forecast_reference_time", 1),
        forecast_period=n_timesteps,
        latitude=z_squeezed.sizes.get("latitude", 1),
        longitude=z_squeezed.sizes.get("longitude", 1),
    )

    return z_expanded


def distribute_daily_precip_to_6hourly(ds_daily, forecast_periods_6h):
    """
    Distribute daily total precipitation to 6-hourly intervals.

    Each daily value is divided by 4 and assigned to the four corresponding
    6-hour intervals.

    Parameters
    ----------
    ds_daily : xarray.Dataset
        Dataset containing daily variables
    forecast_periods_6h : array-like
        The 6-hourly forecast periods (e.g., [6, 12, 18, 24, 30, 36, 42, 48] in hours)

    Returns
    -------
    xarray.DataArray
        tp variable at 6-hourly resolution
    """
    tp_daily = ds_daily["tp"]
    daily_periods_raw = ds_daily["forecast_period"].values

    # Convert to hours
    forecast_periods_6h_hours = forecast_period_to_hours(forecast_periods_6h)
    daily_periods_hours = forecast_period_to_hours(daily_periods_raw)

    # Get dimensions
    n_number = tp_daily.sizes.get("number", tp_daily.shape[0])
    n_ref_time = tp_daily.sizes.get("forecast_reference_time", 1)
    n_lat = tp_daily.sizes.get("latitude", 1)
    n_lon = tp_daily.sizes.get("longitude", 1)
    n_6h = len(forecast_periods_6h)

    # Create output array
    tp_6h = np.zeros((n_number, n_ref_time, n_6h, n_lat, n_lon), dtype=np.float32)

    # For each daily period, distribute to 4 6-hourly periods
    for i, daily_period_hours in enumerate(daily_periods_hours):
        # Find the 4 6-hourly periods that belong to this day
        # Day ending at hour 24 -> 6h intervals at 6, 12, 18, 24
        # Day ending at hour 48 -> 6h intervals at 30, 36, 42, 48
        start_hour = daily_period_hours - 24 + 6  # first 6h period of this day

        # Get indices of the 4 6-hourly periods for this day
        indices_6h = []
        for h in range(4):
            hour = start_hour + h * 6
            idx = np.where(np.isclose(forecast_periods_6h_hours, hour))[0]
            if len(idx) > 0:
                indices_6h.append(idx[0])

        # Divide daily value by 4 and assign to each 6-hourly period
        daily_value = tp_daily.isel(forecast_period=i).values
        for idx in indices_6h:
            tp_6h[:, :, idx, :, :] = daily_value / 4.0

    return tp_6h


def distribute_daily_radiation_to_6hourly_flux(ds_daily, var_name, forecast_periods_6h):
    """
    Distribute daily radiation to 6-hourly intervals and convert to flux.

    For radiation (strd, ssrd):
    - Zero in first 6-hour interval (night: 00-06)
    - Half the daily value in middle two intervals (day: 06-12, 12-18)
    - Zero in last 6-hour interval (night: 18-24)

    The daily accumulated radiation (J/m²) is converted to flux (W/m²) by
    dividing by the time interval (6 hours = 21600 seconds).

    Parameters
    ----------
    ds_daily : xarray.Dataset
        Dataset containing daily variables
    var_name : str
        Variable name ('strd' or 'ssrd')
    forecast_periods_6h : array-like
        The 6-hourly forecast periods

    Returns
    -------
    np.ndarray
        Radiation flux (W/m²) at 6-hourly resolution
    """
    rad_daily = ds_daily[var_name]
    daily_periods_raw = ds_daily["forecast_period"].values

    # Convert to hours
    forecast_periods_6h_hours = forecast_period_to_hours(forecast_periods_6h)
    daily_periods_hours = forecast_period_to_hours(daily_periods_raw)

    # Get dimensions
    n_number = rad_daily.sizes.get("number", rad_daily.shape[0])
    n_ref_time = rad_daily.sizes.get("forecast_reference_time", 1)
    n_lat = rad_daily.sizes.get("latitude", 1)
    n_lon = rad_daily.sizes.get("longitude", 1)
    n_6h = len(forecast_periods_6h)

    # Time interval in seconds (6 hours)
    dt_seconds = 6 * 3600  # 21600 seconds

    # Create output array
    flux_6h = np.zeros((n_number, n_ref_time, n_6h, n_lat, n_lon), dtype=np.float32)

    # For each daily period, distribute to 4 6-hourly periods
    for i, daily_period_hours in enumerate(daily_periods_hours):
        # Find the 4 6-hourly periods that belong to this day
        start_hour = daily_period_hours - 24 + 6

        # Get indices of the 4 6-hourly periods for this day
        hours_6h = [start_hour + h * 6 for h in range(4)]
        indices_6h = []
        for hour in hours_6h:
            idx = np.where(np.isclose(forecast_periods_6h_hours, hour))[0]
            if len(idx) > 0:
                indices_6h.append(idx[0])

        if len(indices_6h) != 4:
            print(
                f"Warning: Expected 4 6-hourly indices for day {daily_period_hours}h, "
                + f"got {len(indices_6h)}"
            )
            continue

        # Daily accumulated value (J/m²)
        daily_value = rad_daily.isel(forecast_period=i).values

        # Convert to flux (W/m²): energy per 6-hour interval / time in seconds
        # Half the daily energy goes to each of the two "day" intervals
        # Flux = (daily_value / 2) / dt_seconds
        flux_value = (daily_value / 2.0) / dt_seconds

        # First interval (night): 0
        flux_6h[:, :, indices_6h[0], :, :] = 0.0
        # Second interval (day): flux
        flux_6h[:, :, indices_6h[1], :, :] = flux_value
        # Third interval (day): flux
        flux_6h[:, :, indices_6h[2], :, :] = flux_value
        # Fourth interval (night): 0
        flux_6h[:, :, indices_6h[3], :, :] = 0.0

    return flux_6h


def main():
    parser = argparse.ArgumentParser(
        description="Convert SEAS5 constant and daily variables to 6-hourly resolution"
    )
    parser.add_argument(
        "--const", required=True, help="Path to netCDF file with constant variables (z)"
    )
    parser.add_argument(
        "--daily",
        required=True,
        help="Path to netCDF file with daily variables (strd, ssrd, tp)",
    )
    parser.add_argument(
        "--hourly", required=True, help="Path to netCDF file with 6-hourly variables"
    )
    parser.add_argument(
        "--output",
        default=None,
        help="Output file path (default: modify hourly file in place)",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Print what would be done without writing output",
    )

    args = parser.parse_args()

    # Determine output file early
    output_file = args.output if args.output else args.hourly
    in_place = output_file == args.hourly

    # Check if output directory exists
    output_dir = os.path.dirname(output_file)
    if output_dir and not os.path.exists(output_dir):
        print(f"Error: Output directory does not exist: {output_dir}")
        print(f"Please create it with: mkdir -p {output_dir}")
        sys.exit(1)

    if in_place:
        print(f"Note: Modifying {args.hourly} in place")

    # Load datasets - use load() to ensure data is in memory
    # This is important for in-place modification
    print(f"Loading constant file: {args.const}")
    ds_const = xr.open_dataset(args.const).load()

    print(f"Loading daily file: {args.daily}")
    ds_daily = xr.open_dataset(args.daily).load()

    print(f"Loading 6-hourly file: {args.hourly}")
    ds_6h = xr.open_dataset(args.hourly).load()

    # Get 6-hourly forecast periods
    forecast_periods_6h = ds_6h["forecast_period"].values
    n_timesteps = len(forecast_periods_6h)

    # Convert to hours for display
    fp_6h_hours = forecast_period_to_hours(forecast_periods_6h)
    fp_daily_hours = forecast_period_to_hours(ds_daily["forecast_period"].values)

    print(f"6-hourly forecast periods (hours): {fp_6h_hours}")
    print(f"Daily forecast periods (hours): {fp_daily_hours}")

    # 1. Expand constant z to 6-hourly
    print("Expanding constant z to 6-hourly resolution...")
    z_const = ds_const["z"]

    # Squeeze out the singleton forecast_period dimension from constant
    z_squeezed = z_const.squeeze("forecast_period", drop=True)

    # Create z_6h by tiling the constant value
    # Shape: (number, forecast_reference_time, forecast_period, latitude, longitude)
    z_data = np.tile(
        z_squeezed.values[:, :, np.newaxis, :, :], (1, 1, n_timesteps, 1, 1)
    )

    # 2. Distribute daily precipitation to 6-hourly
    print("Distributing daily precipitation to 6-hourly...")
    tp_6h = distribute_daily_precip_to_6hourly(ds_daily, forecast_periods_6h)

    # 3. Distribute daily radiation to 6-hourly and convert to flux (W/m²)
    print("Distributing daily strd (thermal radiation) to 6-hourly flux (flds)...")
    flds_6h = distribute_daily_radiation_to_6hourly_flux(
        ds_daily, "strd", forecast_periods_6h
    )

    print("Distributing daily ssrd (solar radiation) to 6-hourly flux (fsds)...")
    fsds_6h = distribute_daily_radiation_to_6hourly_flux(
        ds_daily, "ssrd", forecast_periods_6h
    )

    # Create new dataset with all variables
    print("Creating merged dataset...")

    # Copy the 6-hourly dataset structure
    ds_out = ds_6h.copy(deep=True)

    # Add z variable
    ds_out["z"] = xr.DataArray(
        data=z_data,
        dims=[
            "number",
            "forecast_reference_time",
            "forecast_period",
            "latitude",
            "longitude",
        ],
        attrs=z_const.attrs,
    )

    # Add tp variable
    ds_out["tp"] = xr.DataArray(
        data=tp_6h,
        dims=[
            "number",
            "forecast_reference_time",
            "forecast_period",
            "latitude",
            "longitude",
        ],
        attrs=ds_daily["tp"].attrs,
    )

    # Add flds variable (thermal radiation flux, converted from strd)
    ds_out["flds"] = xr.DataArray(
        data=flds_6h,
        dims=[
            "number",
            "forecast_reference_time",
            "forecast_period",
            "latitude",
            "longitude",
        ],
        attrs={
            "units": "W m-2",
            "long_name": "Downward longwave radiation at surface",
            "standard_name": "surface_downwelling_longwave_flux_in_air",
            "description": "Converted from daily accumulated strd by distributing to 6-hourly intervals and dividing by time",
        },
    )

    # Add fsds variable (solar radiation flux, converted from ssrd)
    ds_out["fsds"] = xr.DataArray(
        data=fsds_6h,
        dims=[
            "number",
            "forecast_reference_time",
            "forecast_period",
            "latitude",
            "longitude",
        ],
        attrs={
            "units": "W m-2",
            "long_name": "Downward shortwave radiation at surface",
            "standard_name": "surface_downwelling_shortwave_flux_in_air",
            "description": "Converted from daily accumulated ssrd by distributing to 6-hourly intervals and dividing by time",
        },
    )

    # Write output (unless dry-run)
    if args.dry_run:
        print(f"\n[DRY RUN] Would write output to: {output_file}")
        print("[DRY RUN] Output dataset structure:")
        print(ds_out)
    else:
        print(f"Writing output to: {output_file}")

        # Build encoding dict preserving original encodings where possible
        encoding = {}

        # Preserve encoding for all coordinates from the 6-hourly file
        for coord in ds_6h.coords:
            if ds_6h[coord].encoding:
                encoding[coord] = {
                    k: v
                    for k, v in ds_6h[coord].encoding.items()
                    if k not in ("source", "original_shape")
                }

        # Preserve encoding for existing data variables from the 6-hourly file
        for var in ds_6h.data_vars:
            if var in ds_out.data_vars and ds_6h[var].encoding:
                encoding[var] = {
                    k: v
                    for k, v in ds_6h[var].encoding.items()
                    if k not in ("source", "original_shape")
                }
                # Ensure compression is enabled
                encoding[var].setdefault("zlib", True)
                encoding[var].setdefault("complevel", 4)

        # Preserve encoding for z from constant file
        if "z" in ds_out.data_vars and ds_const["z"].encoding:
            encoding["z"] = {
                k: v
                for k, v in ds_const["z"].encoding.items()
                if k not in ("source", "original_shape")
            }
            encoding["z"].setdefault("zlib", True)
            encoding["z"].setdefault("complevel", 4)

        # Preserve encoding for tp from daily file
        if "tp" in ds_out.data_vars and ds_daily["tp"].encoding:
            encoding["tp"] = {
                k: v
                for k, v in ds_daily["tp"].encoding.items()
                if k not in ("source", "original_shape")
            }
            encoding["tp"].setdefault("zlib", True)
            encoding["tp"].setdefault("complevel", 4)

        # For flds/fsds (derived from strd/ssrd), use compression with float32
        for var in ("flds", "fsds"):
            if var in ds_out.data_vars and var not in encoding:
                encoding[var] = {"zlib": True, "complevel": 4, "dtype": "float32"}

        ds_out.to_netcdf(output_file, encoding=encoding)

    # Close datasets
    ds_const.close()
    ds_daily.close()
    ds_6h.close()

    print("Done!")

    # Print summary
    print("\nSummary of added variables:")
    print(f"  z: constant orography, broadcast to {n_timesteps} time steps")
    print("  tp: daily precipitation divided by 4 for each 6-hourly interval")
    print(
        "  flds: thermal radiation flux (W/m²) - "
        + "0 at night, flux at day (converted from strd)"
    )
    print(
        "  fsds: solar radiation flux (W/m²) - "
        + "0 at night, flux at day (converted from ssrd)"
    )


if __name__ == "__main__":
    main()
