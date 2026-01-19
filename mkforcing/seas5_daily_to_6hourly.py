#!/usr/bin/env python3
"""
Convert SEAS5 constant, daily, and 6-hourly variables to a target time resolution.

This script:
1. Adds constant `z` (orography) broadcast to all time steps
2. Expands 6-hourly variables (msl, u10, v10, t2m, d2m) to target frequency
3. Converts daily `tp` (total precipitation) by dividing equally across intervals
4. Converts daily `strd` (thermal radiation) to flux, distributed equally
5. Converts daily `ssrd` (solar radiation) to flux with bell-shaped diurnal cycle
   (cosine distribution: zero at 6:00 and 18:00, peak at noon)

Usage:
    python seas5_daily_to_6hourly.py --const <const_file> --daily <daily_file> \\
        --hourly <6h_file> [--output <output_file>] [--frequency <hours>]

If --output is not specified, the output file will be named based on frequency.
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


def generate_target_forecast_periods(input_periods_hours, frequency_hours):
    """
    Generate target forecast periods at the desired frequency.

    Parameters
    ----------
    input_periods_hours : np.ndarray
        Input forecast periods in hours (e.g., [6, 12, 18, 24, 30, 36, 42, 48])
    frequency_hours : int
        Target frequency in hours (e.g., 1 for hourly)

    Returns
    -------
    np.ndarray
        Target forecast periods in hours
    """
    start_hour = frequency_hours  # First period (e.g., 1 for hourly, 6 for 6-hourly)
    end_hour = int(np.max(input_periods_hours))
    return np.arange(start_hour, end_hour + frequency_hours, frequency_hours, dtype=float)


def expand_6hourly_to_target(ds_6h, var_name, target_periods_hours, frequency_hours):
    """
    Expand a 6-hourly variable to target frequency by repeating values.

    Each 6-hourly value is assigned to all sub-intervals within that 6-hour window.
    For example, with hourly output, the value at hour 6 is assigned to hours 1-6.

    Parameters
    ----------
    ds_6h : xarray.Dataset
        Dataset containing 6-hourly variables
    var_name : str
        Variable name to expand
    target_periods_hours : np.ndarray
        Target forecast periods in hours
    frequency_hours : int
        Target frequency in hours

    Returns
    -------
    np.ndarray
        Variable data at target frequency
    """
    var_data = ds_6h[var_name]
    input_periods_hours = forecast_period_to_hours(ds_6h["forecast_period"].values)

    # Get dimensions
    n_number = var_data.sizes.get("number", var_data.shape[0])
    n_ref_time = var_data.sizes.get("forecast_reference_time", 1)
    n_lat = var_data.sizes.get("latitude", 1)
    n_lon = var_data.sizes.get("longitude", 1)
    n_target = len(target_periods_hours)

    # Create output array
    output = np.zeros((n_number, n_ref_time, n_target, n_lat, n_lon), dtype=np.float32)

    # For each target period, find the corresponding 6-hourly period
    for i, target_hour in enumerate(target_periods_hours):
        # Find the 6-hourly period that contains this target hour
        # Hour 1-6 -> 6h period at hour 6, Hour 7-12 -> 6h period at hour 12, etc.
        containing_6h = int(np.ceil(target_hour / 6.0) * 6)
        idx_6h = np.where(np.isclose(input_periods_hours, containing_6h))[0]
        if len(idx_6h) > 0:
            output[:, :, i, :, :] = var_data.isel(forecast_period=idx_6h[0]).values

    return output


def distribute_daily_to_target(ds_daily, var_name, target_periods_hours, frequency_hours):
    """
    Distribute daily accumulated variable to target frequency intervals.

    Each daily value is divided equally among all intervals within that day.

    Parameters
    ----------
    ds_daily : xarray.Dataset
        Dataset containing daily variables
    var_name : str
        Variable name (e.g., 'tp')
    target_periods_hours : np.ndarray
        Target forecast periods in hours
    frequency_hours : int
        Target frequency in hours

    Returns
    -------
    np.ndarray
        Variable data at target frequency
    """
    var_daily = ds_daily[var_name]
    daily_periods_hours = forecast_period_to_hours(ds_daily["forecast_period"].values)

    # Get dimensions
    n_number = var_daily.sizes.get("number", var_daily.shape[0])
    n_ref_time = var_daily.sizes.get("forecast_reference_time", 1)
    n_lat = var_daily.sizes.get("latitude", 1)
    n_lon = var_daily.sizes.get("longitude", 1)
    n_target = len(target_periods_hours)

    # Number of intervals per day
    intervals_per_day = int(24 / frequency_hours)

    # Create output array
    output = np.zeros((n_number, n_ref_time, n_target, n_lat, n_lon), dtype=np.float32)

    # For each daily period, distribute to target intervals
    for i, daily_period_hours in enumerate(daily_periods_hours):
        # Find target intervals belonging to this day
        day_start = daily_period_hours - 24 + frequency_hours
        day_end = daily_period_hours

        indices = []
        for j, target_hour in enumerate(target_periods_hours):
            if day_start <= target_hour <= day_end:
                indices.append(j)

        if len(indices) > 0:
            # Divide daily value by number of intervals
            daily_value = var_daily.isel(forecast_period=i).values
            for idx in indices:
                output[:, :, idx, :, :] = daily_value / len(indices)

    return output


def distribute_thermal_radiation_to_flux(ds_daily, var_name, target_periods_hours, frequency_hours):
    """
    Distribute daily thermal radiation to target frequency and convert to flux.

    Thermal radiation is distributed evenly across all intervals (day and night).

    Parameters
    ----------
    ds_daily : xarray.Dataset
        Dataset containing daily variables
    var_name : str
        Variable name (typically 'strd')
    target_periods_hours : np.ndarray
        Target forecast periods in hours
    frequency_hours : int
        Target frequency in hours

    Returns
    -------
    np.ndarray
        Radiation flux (W/m²) at target frequency
    """
    rad_daily = ds_daily[var_name]
    daily_periods_hours = forecast_period_to_hours(ds_daily["forecast_period"].values)

    # Get dimensions
    n_number = rad_daily.sizes.get("number", rad_daily.shape[0])
    n_ref_time = rad_daily.sizes.get("forecast_reference_time", 1)
    n_lat = rad_daily.sizes.get("latitude", 1)
    n_lon = rad_daily.sizes.get("longitude", 1)
    n_target = len(target_periods_hours)

    # Time interval in seconds
    dt_seconds = frequency_hours * 3600

    # Number of intervals per day
    intervals_per_day = int(24 / frequency_hours)

    # Create output array
    flux = np.zeros((n_number, n_ref_time, n_target, n_lat, n_lon), dtype=np.float32)

    # For each daily period, distribute to target intervals
    for i, daily_period_hours in enumerate(daily_periods_hours):
        day_start = daily_period_hours - 24 + frequency_hours
        day_end = daily_period_hours

        indices = []
        for j, target_hour in enumerate(target_periods_hours):
            if day_start <= target_hour <= day_end:
                indices.append(j)

        if len(indices) > 0:
            # Daily accumulated value (J/m²)
            daily_value = rad_daily.isel(forecast_period=i).values
            # Flux = (daily_value / n_intervals) / dt_seconds
            flux_value = (daily_value / len(indices)) / dt_seconds
            for idx in indices:
                flux[:, :, idx, :, :] = flux_value

    return flux


def distribute_solar_radiation_to_flux(ds_daily, var_name, target_periods_hours, frequency_hours):
    """
    Distribute daily solar radiation with bell-shaped diurnal cycle.

    Uses a cosine distribution:
    - Zero radiation between 18:00 and 06:00 (night)
    - Bell-shaped curve between 06:00 and 18:00, peak at noon
    - weight(h) = max(0, cos((h - 12) * π / 12))

    Parameters
    ----------
    ds_daily : xarray.Dataset
        Dataset containing daily variables
    var_name : str
        Variable name (typically 'ssrd')
    target_periods_hours : np.ndarray
        Target forecast periods in hours
    frequency_hours : int
        Target frequency in hours

    Returns
    -------
    np.ndarray
        Radiation flux (W/m²) at target frequency
    """
    rad_daily = ds_daily[var_name]
    daily_periods_hours = forecast_period_to_hours(ds_daily["forecast_period"].values)

    # Get dimensions
    n_number = rad_daily.sizes.get("number", rad_daily.shape[0])
    n_ref_time = rad_daily.sizes.get("forecast_reference_time", 1)
    n_lat = rad_daily.sizes.get("latitude", 1)
    n_lon = rad_daily.sizes.get("longitude", 1)
    n_target = len(target_periods_hours)

    # Time interval in seconds
    dt_seconds = frequency_hours * 3600

    # Create output array
    flux = np.zeros((n_number, n_ref_time, n_target, n_lat, n_lon), dtype=np.float32)

    # For each daily period, distribute with bell-shaped weights
    for i, daily_period_hours in enumerate(daily_periods_hours):
        day_start = daily_period_hours - 24 + frequency_hours
        day_end = daily_period_hours

        # Find indices and calculate weights for this day
        indices = []
        weights = []
        for j, target_hour in enumerate(target_periods_hours):
            if day_start <= target_hour <= day_end:
                indices.append(j)
                # Hour of day (0-24) - use center of interval
                hour_of_day = ((target_hour - frequency_hours / 2) % 24)
                # Bell-shaped weight: cosine centered at noon
                # cos((h - 12) * π / 12) gives 0 at h=6 and h=18, 1 at h=12
                weight = max(0.0, np.cos((hour_of_day - 12) * np.pi / 12))
                weights.append(weight)

        if len(indices) > 0 and sum(weights) > 0:
            # Normalize weights to sum to 1
            weights = np.array(weights)
            weights = weights / weights.sum()

            # Daily accumulated value (J/m²)
            daily_value = rad_daily.isel(forecast_period=i).values

            # Distribute according to weights and convert to flux
            for idx, weight in zip(indices, weights):
                # Energy for this interval = daily_value * weight
                # Flux = energy / dt_seconds
                flux[:, :, idx, :, :] = (daily_value * weight) / dt_seconds

    return flux


def main():
    parser = argparse.ArgumentParser(
        description="Convert SEAS5 variables to target time resolution"
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
        help="Output file path (default: auto-generated based on frequency)",
    )
    parser.add_argument(
        "--frequency",
        type=int,
        default=6,
        choices=[1, 2, 3, 6],
        help="Target frequency in hours (default: 6). Must divide 24 evenly.",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Print what would be done without writing output",
    )

    args = parser.parse_args()

    frequency_hours = args.frequency

    # Determine output file
    if args.output:
        output_file = args.output
    else:
        base = os.path.splitext(args.hourly)[0]
        output_file = f"{base}_{frequency_hours}h.nc"

    # Check if output directory exists
    output_dir = os.path.dirname(output_file)
    if output_dir and not os.path.exists(output_dir):
        print(f"Error: Output directory does not exist: {output_dir}")
        print(f"Please create it with: mkdir -p {output_dir}")
        sys.exit(1)

    print(f"Target frequency: {frequency_hours} hours")
    print(f"Output file: {output_file}")

    # Load datasets
    print(f"Loading constant file: {args.const}")
    ds_const = xr.open_dataset(args.const).load()

    print(f"Loading daily file: {args.daily}")
    ds_daily = xr.open_dataset(args.daily).load()

    print(f"Loading 6-hourly file: {args.hourly}")
    ds_6h = xr.open_dataset(args.hourly).load()

    # Get input periods and generate target periods
    input_periods_hours = forecast_period_to_hours(ds_6h["forecast_period"].values)
    target_periods_hours = generate_target_forecast_periods(input_periods_hours, frequency_hours)
    n_target = len(target_periods_hours)

    daily_periods_hours = forecast_period_to_hours(ds_daily["forecast_period"].values)

    print(f"Input 6-hourly periods (hours): {input_periods_hours}")
    print(f"Daily periods (hours): {daily_periods_hours}")
    print(f"Target periods (hours): {target_periods_hours}")

    # 1. Expand constant z
    print("Expanding constant z...")
    z_const = ds_const["z"]
    z_squeezed = z_const.squeeze("forecast_period", drop=True)
    z_data = np.tile(
        z_squeezed.values[:, :, np.newaxis, :, :], (1, 1, n_target, 1, 1)
    )

    # 2. Expand 6-hourly variables
    vars_6h = ["msl", "u10", "v10", "t2m", "d2m"]
    expanded_vars = {}
    for var in vars_6h:
        if var in ds_6h.data_vars:
            print(f"Expanding {var} to {frequency_hours}-hourly...")
            expanded_vars[var] = expand_6hourly_to_target(
                ds_6h, var, target_periods_hours, frequency_hours
            )

    # 3. Distribute daily precipitation
    print("Distributing daily precipitation...")
    tp_target = distribute_daily_to_target(
        ds_daily, "tp", target_periods_hours, frequency_hours
    )

    # 4. Distribute thermal radiation (evenly)
    print("Distributing thermal radiation (flds)...")
    flds_target = distribute_thermal_radiation_to_flux(
        ds_daily, "strd", target_periods_hours, frequency_hours
    )

    # 5. Distribute solar radiation (bell-shaped)
    print("Distributing solar radiation with bell-shaped diurnal cycle (fsds)...")
    fsds_target = distribute_solar_radiation_to_flux(
        ds_daily, "ssrd", target_periods_hours, frequency_hours
    )

    # Create output dataset
    print("Creating output dataset...")

    dims = ["number", "forecast_reference_time", "forecast_period", "latitude", "longitude"]

    # Start with coordinates
    ds_out = xr.Dataset(
        coords={
            "number": ds_6h["number"],
            "forecast_reference_time": ds_6h["forecast_reference_time"],
            "forecast_period": target_periods_hours,
            "latitude": ds_6h["latitude"],
            "longitude": ds_6h["longitude"],
        }
    )

    # Copy coordinate attributes
    ds_out["forecast_period"].attrs = {
        "long_name": "time since forecast_reference_time",
        "standard_name": "forecast_period",
        "units": "hours",
    }

    # Compute valid_time = forecast_reference_time + forecast_period
    # forecast_reference_time may be datetime64 or seconds since epoch
    ref_time = ds_6h["forecast_reference_time"].values[0]
    if np.issubdtype(type(ref_time), np.datetime64):
        # Convert datetime64 to seconds since epoch
        ref_time_seconds = (ref_time - np.datetime64('1970-01-01T00:00:00')) / np.timedelta64(1, 's')
    else:
        ref_time_seconds = float(ref_time)
    valid_time_seconds = (ref_time_seconds + target_periods_hours * 3600).astype(np.int64)
    # ds_out["valid_time"] = xr.DataArray(
    #     data=valid_time_seconds,
    #     dims=["forecast_period"],
    #     attrs={
    #         "standard_name": "time",
    #         "long_name": "time",
    #         "units": "seconds since 1970-01-01",
    #         "calendar": "proleptic_gregorian",
    #     },
    # )

    # Also add time in hours since 1900-01-01 (common format for climate data)
    # Hours from 1900-01-01 to 1970-01-01: 613608 hours
    hours_1900_to_1970 = 613608
    time_hours = (valid_time_seconds / 3600 + hours_1900_to_1970).astype(np.int32)
    ds_out["valid_time"] = xr.DataArray(
        data=time_hours,
        dims=["forecast_period"],
        attrs={
            "standard_name": "time",
            "long_name": "time",
            "units": "hours since 1900-01-01 00:00:00.0",
            "calendar": "gregorian",
            "axis": "T",
        },
    )

    # Add z
    ds_out["z"] = xr.DataArray(data=z_data, dims=dims, attrs=z_const.attrs)

    # Add expanded 6-hourly variables
    for var, data in expanded_vars.items():
        ds_out[var] = xr.DataArray(data=data, dims=dims, attrs=ds_6h[var].attrs)

    # Add precipitation
    ds_out["tp"] = xr.DataArray(data=tp_target, dims=dims, attrs=ds_daily["tp"].attrs)

    # Add radiation fluxes
    ds_out["flds"] = xr.DataArray(
        data=flds_target,
        dims=dims,
        attrs={
            "units": "W m-2",
            "long_name": "Downward longwave radiation at surface",
            "standard_name": "surface_downwelling_longwave_flux_in_air",
            "description": f"Converted from daily strd, distributed equally to {frequency_hours}-hourly intervals",
        },
    )

    ds_out["fsds"] = xr.DataArray(
        data=fsds_target,
        dims=dims,
        attrs={
            "units": "W m-2",
            "long_name": "Downward shortwave radiation at surface",
            "standard_name": "surface_downwelling_shortwave_flux_in_air",
            "description": f"Converted from daily ssrd with bell-shaped diurnal cycle (cosine, peak at noon)",
        },
    )

    # Copy global attributes
    ds_out.attrs = ds_6h.attrs.copy()
    ds_out.attrs["frequency"] = f"{frequency_hours} hours"

    # Write output
    if args.dry_run:
        print(f"\n[DRY RUN] Would write output to: {output_file}")
        print("[DRY RUN] Output dataset structure:")
        print(ds_out)
    else:
        print(f"Writing output to: {output_file}")

        # Build encoding
        valid_encodings = {
            "compression", "dtype", "least_significant_digit", "zlib",
            "_FillValue", "fletcher32", "complevel", "chunksizes",
            "shuffle", "contiguous", "calendar", "units",
        }

        def filter_encoding(enc):
            return {k: v for k, v in enc.items() if k in valid_encodings}

        encoding = {}

        # Default encoding for all data variables
        for var in ds_out.data_vars:
            encoding[var] = {"zlib": True, "complevel": 4, "dtype": "float32"}

        # Preserve encoding from source files where applicable
        for var in vars_6h:
            if var in ds_6h.data_vars and ds_6h[var].encoding:
                encoding[var].update(filter_encoding(ds_6h[var].encoding))

        if ds_const["z"].encoding:
            encoding["z"].update(filter_encoding(ds_const["z"].encoding))

        if ds_daily["tp"].encoding:
            encoding["tp"].update(filter_encoding(ds_daily["tp"].encoding))

        # valid_time should be int64, time should be int32
        # encoding["valid_time"] = {"dtype": "int64"}
        encoding["valid_time"] = {"dtype": "int32"}

        ds_out.to_netcdf(output_file, encoding=encoding)

    # Close datasets
    ds_const.close()
    ds_daily.close()
    ds_6h.close()

    print("Done!")

    # Print summary
    print(f"\nSummary (target frequency: {frequency_hours} hours):")
    print(f"  z: constant orography, broadcast to {n_target} time steps")
    for var in vars_6h:
        if var in expanded_vars:
            print(f"  {var}: expanded from 6-hourly by value repetition")
    print(f"  tp: daily precipitation distributed equally")
    print(f"  flds: thermal radiation flux - distributed equally (day and night)")
    print(f"  fsds: solar radiation flux - bell-shaped (cosine, 6-18h, peak at noon)")


if __name__ == "__main__":
    main()
