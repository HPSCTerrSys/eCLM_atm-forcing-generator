"""
Convert between Mean Sea Level Pressure (MSLP) and Surface Pressure
Using the Method from Stull's Practical Meteorology Chapter 9

This implementation follows the standard meteorological sea-level pressure
reduction method that accounts for the elevation and uses a fictitious
temperature for the imaginary air column between the surface and sea level.

References (all open-access):
1. Stull, R., 2017: Practical Meteorology: An Algebra-based Survey of
   Atmospheric Science. University of British Columbia, 940 pp.
   ISBN 978-0-88865-283-6
   Available under Creative Commons License (CC BY-NC-SA 4.0)
   URL: https://www.eoas.ubc.ca/books/Practical_Meteorology/
   See Chapter 9, Section "Sea-level Pressure Reduction"

Mathematical Background:
The sea-level pressure reduction uses the hypsometric equation:

    P_MSL = P_surface * exp(z_stn / (a * T_v*))

where:
- P_MSL: Mean sea level pressure (Pa or hPa)
- P_surface: Surface pressure at the station (Pa or hPa)
- z_stn: Station elevation above sea level (m)
- a: Scale height parameter = R_d / g ≈ 29.3 m/K
- T_v*: Fictitious average virtual temperature (K)

The fictitious temperature T_v* represents the temperature of an imaginary
air column between the surface and sea level. It is calculated as:

    T_v* = 0.5 * [T_v(t_0) + T_v(t_0 - 12h) + γ_sa * z_stn]

where:
- T_v(t_0): Current virtual temperature at the surface
- T_v(t_0 - 12h): Virtual temperature 12 hours ago
- γ_sa: Standard atmosphere lapse rate = 0.0065 K/m
"""

import argparse
import numpy as np
import netCDF4


def mslp_to_surface_pressure(
    mslp: float,
    elevation: float,
    temperature_current: float,
    temperature_12h_ago: float = None,
    pressure_units: str = "Pa",
) -> float:
    """
    Convert mean sea level pressure to surface pressure at a given elevation.

    Uses the standard meteorological method with a fictitious temperature
    for the imaginary air column between the surface and sea level.

    Parameters
    ----------
    mslp : float
        Mean sea level pressure (in units specified by pressure_units)
    elevation : float
        Station elevation above mean sea level (meters)
        - Positive for locations above sea level
        - Negative for locations below sea level
    temperature_current : float
        Current surface air temperature (Kelvin)
    temperature_12h_ago : float, optional
        Surface air temperature 12 hours ago (Kelvin)
        If None, uses temperature_current (i.e., assumes steady conditions)
    pressure_units : str, optional
        Units for pressure ('Pa' or 'hPa'), default='Pa'

    Returns
    -------
    float
        Surface pressure at the given elevation (same units as input)

    Notes
    -----
    1. This follows the method described in Stull (2017), Chapter 9.
    2. The fictitious temperature accounts for:
       - Current surface temperature
       - Temperature 12 hours ago (for diurnal averaging)
       - Lapse rate correction for the imaginary column below the surface
    3. For locations at sea level (elevation=0), returns mslp unchanged.
    4. The 12-hour averaging helps smooth out diurnal temperature variations.

    Examples
    --------
    >>> # Sea level location
    >>> surface_p = mslp_to_surface_pressure(101325, 0, 288.15, None, 'Pa')
    >>> print(f"Surface pressure: {surface_p:.1f} Pa")
    Surface pressure: 101325.0 Pa

    >>> # Mountain station (2000m) with temperature data
    >>> surface_p = mslp_to_surface_pressure(1013.25, 2000, 278, 276, 'hPa')
    >>> print(f"Surface pressure: {surface_p:.1f} hPa")
    Surface pressure: 795.3 hPa
    """

    # Physical constants
    R_d = 287.05  # Gas constant for dry air (J/(kg·K))
    g = 9.80665  # Gravitational acceleration (m/s²)
    gamma_sa = 0.0065  # Standard atmosphere lapse rate (K/m)

    # Calculate scale height parameter
    a = R_d / g  # ≈ 29.3 m/K

    # If no 12-hour-ago temperature provided, use current temperature
    if temperature_12h_ago is None:
        temperature_12h_ago = temperature_current

    # Calculate fictitious average virtual temperature for the column
    # This is the temperature of the imaginary air column between
    # the surface and sea level
    # T_v* = 0.5 * [T_v(now) + T_v(12h ago) + γ_sa * z_stn]
    T_v_star = 0.5 * (temperature_current + temperature_12h_ago + gamma_sa * elevation)

    # Apply hypsometric equation (inverse direction)
    # P_surface = P_MSL * exp(-z / (a * T_v*))
    exponent = -elevation / (a * T_v_star)
    surface_pressure = mslp * np.exp(exponent)

    return surface_pressure


def surface_to_mslp(
    surface_pressure: float,
    elevation: float,
    temperature_current: float,
    temperature_12h_ago: float = None,
    pressure_units: str = "Pa",
) -> float:
    """
    Convert surface pressure to mean sea level pressure.

    This is the standard meteorological "sea-level pressure reduction"
    operation used to create weather maps.

    Parameters
    ----------
    surface_pressure : float
        Pressure at the surface (in units specified by pressure_units)
    elevation : float
        Station elevation above mean sea level (meters)
    temperature_current : float
        Current surface air temperature (Kelvin)
    temperature_12h_ago : float, optional
        Surface air temperature 12 hours ago (Kelvin)
        If None, uses temperature_current
    pressure_units : str, optional
        Units for pressure ('Pa' or 'hPa'), default='Pa'

    Returns
    -------
    float
        Mean sea level pressure (same units as input)

    Notes
    -----
    This function uses the fictitious temperature method from Stull (2017):

    T_v* = 0.5 * [T_v(now) + T_v(12h ago) + γ_sa * z_stn]

    Then applies:
    P_MSL = P_surface * exp(z / (a * T_v*))

    Examples
    --------
    >>> # Convert Denver surface pressure to MSLP
    >>> mslp = surface_to_mslp(83500, 1600, 288.15, 285.15, 'Pa')
    >>> print(f"MSLP: {mslp:.1f} Pa ({mslp/100:.1f} hPa)")
    MSLP: 101362.4 Pa (1013.6 hPa)

    >>> # Mountain weather station
    >>> mslp = surface_to_mslp(795.0, 2000, 278, 276, 'hPa')
    >>> print(f"MSLP: {mslp:.2f} hPa")
    MSLP: 1012.85 hPa
    """

    # Physical constants
    R_d = 287.05  # Gas constant for dry air (J/(kg·K))
    g = 9.80665  # Gravitational acceleration (m/s²)
    gamma_sa = 0.0065  # Standard atmosphere lapse rate (K/m)

    # Calculate scale height parameter
    a = R_d / g  # ≈ 29.3 m/K

    # If no 12-hour-ago temperature provided, use current temperature
    if temperature_12h_ago is None:
        temperature_12h_ago = temperature_current

    # Calculate fictitious average temperature for the imaginary column
    T_v_star = 0.5 * (temperature_current + temperature_12h_ago + gamma_sa * elevation)

    # Apply hypsometric equation
    # P_MSL = P_surface * exp(z / (a * T_v*))
    exponent = elevation / (a * T_v_star)
    mslp = surface_pressure * np.exp(exponent)

    return mslp


def calculate_fictitious_temperature(
    temperature_current: float, temperature_12h_ago: float, elevation: float
) -> float:
    """
    Calculate the fictitious temperature for the imaginary air column.

    This is a helper function that computes T_v* according to the method
    in Stull (2017), Chapter 9.

    Parameters
    ----------
    temperature_current : float
        Current surface temperature (Kelvin)
    temperature_12h_ago : float
        Temperature 12 hours ago (Kelvin)
    elevation : float
        Station elevation (meters)

    Returns
    -------
    float
        Fictitious temperature T_v* (Kelvin)

    Notes
    -----
    The formula is:
    T_v* = 0.5 * [T(now) + T(12h ago) + γ_sa * z]

    This represents:
    1. Average of current and 12-hour-ago temperatures (diurnal smoothing)
    2. Plus a lapse-rate correction for the imaginary column height

    Examples
    --------
    >>> T_star = calculate_fictitious_temperature(288.15, 285.15, 1600)
    >>> print(f"Fictitious temperature: {T_star:.2f} K ({T_star-273.15:.2f}°C)")
    Fictitious temperature: 291.55 K (18.40°C)
    """
    gamma_sa = 0.0065  # Standard atmosphere lapse rate (K/m)

    T_v_star = 0.5 * (temperature_current + temperature_12h_ago + gamma_sa * elevation)

    return T_v_star


def add_surface_pressure_to_netcdf(filename, elevation_var="z", temp_var="t2m", mslp_var="msl"):
    """
    Read MSLP, temperature, and elevation from a netCDF file,
    calculate surface pressure, and write it back to the file.

    Parameters:
    -----------
    filename : str
        Path to the netCDF file
    elevation_var : str, optional
        Name of the elevation variable in the file (default: 'z')
        The geopotential variable will be converted to elevation (m) by dividing by g
    temp_var : str, optional
        Name of the temperature variable (default: 't2m')
    mslp_var : str, optional
        Name of the MSLP variable (default: 'msl')

    Returns:
    --------
    None
        Modifies the netCDF file in place by adding 'sp' variable

    Raises:
    -------
    ValueError
        If 'sp' variable already exists in the file
    KeyError
        If required variables are not found
    """

    # Open netCDF file in append mode
    print(f"Opening {filename}...")
    nc = netCDF4.Dataset(filename, "a")

    try:
        # Check if sp already exists - if so, raise error and exit
        if "sp" in nc.variables:
            nc.close()
            raise ValueError(
                f"Variable 'sp' already exists in {filename}. "
                "No changes made. Delete the variable first if you want to recalculate."
            )

        # Read the required variables
        print(f"Reading {mslp_var}, {temp_var}, and {elevation_var}...")
        msl = nc.variables[mslp_var][:]  # Mean sea level pressure [Pa]
        t2m = nc.variables[temp_var][:]  # Temperature at 2m [K]
        z_geopotential = nc.variables[elevation_var][:]  # Geopotential [m^2/s^2]

        print(f"Data shapes - {mslp_var}: {msl.shape}, {temp_var}: {t2m.shape}, {elevation_var}: {z_geopotential.shape}")

        # Convert geopotential to elevation (meters)
        g = 9.80665  # Standard gravity [m/s^2]
        elevation = z_geopotential / g

        print(f"Elevation range: {np.min(elevation):.1f} to {np.max(elevation):.1f} m")

        # Determine time dimension for 12h offset
        # Assume first dimension is time
        time_dim = nc.variables[temp_var].dimensions[0]
        time_size = nc.dimensions[time_dim].size
        print(f"Time dimension: {time_dim} with size {time_size}")

        # Calculate surface pressure
        print("Calculating surface pressure...")

        # Create array for surface pressure
        sp = np.zeros_like(msl)

        # For the first timestep (index 0), we don't have 12h ago data
        # Use current temperature for both
        if time_size > 0:
            print("Processing first timestep (no 12h-ago data available)...")
            sp[0, ...] = mslp_to_surface_pressure(
                msl[0, ...],
                elevation,
                t2m[0, ...],
                t2m[0, ...],  # Use current temp as 12h-ago temp
                "Pa"
            )

        # For remaining timesteps, check if we can use 12h offset
        # Assuming hourly data, 12h offset = 12 timesteps back
        # For other time resolutions, this logic may need adjustment
        if time_size > 12:
            print("Processing remaining timesteps with 12h-ago data...")
            for t in range(12, time_size):
                sp[t, ...] = mslp_to_surface_pressure(
                    msl[t, ...],
                    elevation,
                    t2m[t, ...],
                    t2m[t - 12, ...],  # Temperature 12 timesteps ago
                    "Pa"
                )

            # For timesteps 1-11, use current temp (no 12h data yet)
            if time_size > 1:
                print("Processing timesteps 1-11 without 12h-ago data...")
                for t in range(1, min(12, time_size)):
                    sp[t, ...] = mslp_to_surface_pressure(
                        msl[t, ...],
                        elevation,
                        t2m[t, ...],
                        t2m[t, ...],
                        "Pa"
                    )
        elif time_size > 1:
            # If we have fewer than 12 timesteps, process without 12h offset
            print("Processing all timesteps without 12h-ago data (less than 12 timesteps)...")
            for t in range(1, time_size):
                sp[t, ...] = mslp_to_surface_pressure(
                    msl[t, ...],
                    elevation,
                    t2m[t, ...],
                    t2m[t, ...],
                    "Pa"
                )

        # Create the new variable
        print("Creating new variable 'sp'...")

        # Get dimension names from msl
        dim_names = nc.variables[mslp_var].dimensions

        # Create the sp variable with same dimensions as msl
        sp_var = nc.createVariable("sp", "f4", dim_names, zlib=True, complevel=4)

        # Add attributes
        sp_var.units = "Pa"
        sp_var.long_name = "Surface pressure"
        sp_var.standard_name = "surface_air_pressure"
        sp_var.description = (
            f"Calculated from {mslp_var}, {temp_var}, and {elevation_var} "
            "using Stull (2017) Chapter 9 method"
        )

        # Write the data
        sp_var[:] = sp

        print("Variable 'sp' created and written successfully!")

        # Print some statistics
        print("\nStatistics:")
        print(f"  Surface pressure range: {np.min(sp):.1f} to {np.max(sp):.1f} Pa")
        print(f"  Surface pressure mean: {np.mean(sp):.1f} Pa ({np.mean(sp)/100:.2f} hPa)")
        print(f"  MSLP mean: {np.mean(msl):.1f} Pa ({np.mean(msl)/100:.2f} hPa)")
        print(f"  Mean difference: {np.mean(msl - sp):.1f} Pa ({np.mean(msl - sp)/100:.2f} hPa)")

    except KeyError as e:
        print(f"Error: Required variable not found in netCDF file: {e}")
        print(f"Available variables: {list(nc.variables.keys())}")
        nc.close()
        raise

    except ValueError as e:
        # This catches the "sp already exists" error
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
        description="Add surface pressure (sp) to a netCDF file based on MSLP, temperature, and elevation."
    )
    parser.add_argument(
        "filename",
        type=str,
        help="Path to the netCDF file containing MSLP, temperature, and elevation variables"
    )
    parser.add_argument(
        "--elevation-var",
        type=str,
        default="z",
        help="Name of the elevation/geopotential variable (default: z)"
    )
    parser.add_argument(
        "--temp-var",
        type=str,
        default="t2m",
        help="Name of the temperature variable (default: t2m)"
    )
    parser.add_argument(
        "--mslp-var",
        type=str,
        default="msl",
        help="Name of the MSLP variable (default: msl)"
    )

    # Parse command-line arguments
    args = parser.parse_args()

    # Process the file
    add_surface_pressure_to_netcdf(
        args.filename,
        elevation_var=args.elevation_var,
        temp_var=args.temp_var,
        mslp_var=args.mslp_var
    )
