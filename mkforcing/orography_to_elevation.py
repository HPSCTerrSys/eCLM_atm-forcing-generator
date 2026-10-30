"""
Convert geopotential (orography) to elevation above mean sea level.

Geopotential represents the gravitational potential energy per unit mass
at a given height. It is related to elevation through the gravitational
acceleration constant.

Mathematical Background:
The relationship between geopotential and geometric height is:

    z = Φ / g

where:
- z: Geometric height (elevation) above mean sea level (m)
- Φ: Geopotential (m²/s²)
- g: Standard gravitational acceleration = 9.80665 m/s²

This conversion assumes a constant gravitational acceleration, which is
appropriate for typical atmospheric applications at Earth's surface.

References:
1. WMO Guide to Meteorological Instruments and Methods of Observation
   (WMO-No. 8), 2018 edition
2. ECMWF IFS Documentation - Part III: Dynamics and Numerical Procedures
"""

import argparse
import numpy as np
import netCDF4


def geopotential_to_elevation(geopotential, g=9.80665):
    """
    Convert geopotential to elevation above mean sea level.

    Parameters
    ----------
    geopotential : float or array
        Geopotential [m²/s²]
    g : float, optional
        Standard gravitational acceleration [m/s²] (default: 9.80665)

    Returns
    -------
    elevation : float or array
        Elevation above mean sea level [m]

    Notes
    -----
    The standard gravitational acceleration g = 9.80665 m/s² is defined
    by ISO 80000-3:2006 and is used in meteorological applications.

    Examples
    --------
    >>> # Sea level
    >>> elevation = geopotential_to_elevation(0.0)
    >>> print(f"Elevation: {elevation:.1f} m")
    Elevation: 0.0 m

    >>> # Typical mountain height
    >>> geopotential = 19613.3  # m²/s²
    >>> elevation = geopotential_to_elevation(geopotential)
    >>> print(f"Elevation: {elevation:.1f} m")
    Elevation: 2000.0 m
    """
    elevation = geopotential / g
    return elevation


def elevation_to_geopotential(elevation, g=9.80665):
    """
    Convert elevation to geopotential.

    Parameters
    ----------
    elevation : float or array
        Elevation above mean sea level [m]
    g : float, optional
        Standard gravitational acceleration [m/s²] (default: 9.80665)

    Returns
    -------
    geopotential : float or array
        Geopotential [m²/s²]

    Notes
    -----
    This is the inverse operation of geopotential_to_elevation.

    Examples
    --------
    >>> # Mount Everest
    >>> geopotential = elevation_to_geopotential(8849)
    >>> print(f"Geopotential: {geopotential:.1f} m²/s²")
    Geopotential: 86763.3 m²/s²

    >>> # Below sea level (Dead Sea)
    >>> geopotential = elevation_to_geopotential(-430)
    >>> print(f"Geopotential: {geopotential:.1f} m²/s²")
    Geopotential: -4216.6 m²/s²
    """
    geopotential = elevation * g
    return geopotential


def add_elevation_to_netcdf(filename, geopotential_var="z", elevation_var_name="elevation"):
    """
    Read geopotential (orography) from a netCDF file,
    calculate elevation, and write it back to the file.

    Parameters:
    -----------
    filename : str
        Path to the netCDF file
    geopotential_var : str, optional
        Name of the geopotential variable in the file (default: 'z')
    elevation_var_name : str, optional
        Name for the output elevation variable (default: 'elevation')

    Returns:
    --------
    None
        Modifies the netCDF file in place by adding elevation variable

    Raises:
    -------
    ValueError
        If elevation variable already exists in the file
    KeyError
        If required geopotential variable is not found
    """

    # Open netCDF file in append mode
    print(f"Opening {filename}...")
    nc = netCDF4.Dataset(filename, "a")

    try:
        # Check if elevation variable already exists - if so, raise error and exit
        if elevation_var_name in nc.variables:
            nc.close()
            raise ValueError(
                f"Variable '{elevation_var_name}' already exists in {filename}. "
                "No changes made. Delete the variable first if you want to recalculate."
            )

        # Read the geopotential variable
        print(f"Reading geopotential variable '{geopotential_var}'...")
        geopotential = nc.variables[geopotential_var][:]  # Geopotential [m²/s²]

        print(f"Data shape - {geopotential_var}: {geopotential.shape}")
        print(f"Geopotential range: {np.min(geopotential):.1f} to {np.max(geopotential):.1f} m²/s²")

        # Calculate elevation
        print("Calculating elevation...")
        elevation = geopotential_to_elevation(geopotential)

        print(f"Elevation range: {np.min(elevation):.1f} to {np.max(elevation):.1f} m")

        # Create the new variable
        print(f"Creating new variable '{elevation_var_name}'...")

        # Get dimension names from geopotential variable
        dim_names = nc.variables[geopotential_var].dimensions

        # Create the elevation variable with same dimensions as geopotential
        elevation_var = nc.createVariable(
            elevation_var_name, "f4", dim_names, zlib=True, complevel=4
        )

        # Add attributes
        elevation_var.units = "m"
        elevation_var.long_name = "Elevation above mean sea level"
        elevation_var.standard_name = "surface_altitude"
        elevation_var.description = (
            f"Calculated from geopotential variable '{geopotential_var}' "
            "using z = Φ / g with g = 9.80665 m/s²"
        )

        # Write the data
        elevation_var[:] = elevation

        print(f"Variable '{elevation_var_name}' created and written successfully!")

        # Print some statistics
        print("\nStatistics:")
        print(f"  Elevation range: {np.min(elevation):.1f} to {np.max(elevation):.1f} m")
        print(f"  Elevation mean: {np.mean(elevation):.1f} m")
        print(f"  Elevation std: {np.std(elevation):.1f} m")

        # Identify any interesting features
        if np.min(elevation) < 0:
            print(f"  Below sea level: Yes (min: {np.min(elevation):.1f} m)")
        if np.max(elevation) > 1000:
            print(f"  High elevation areas: Yes (max: {np.max(elevation):.1f} m)")

    except KeyError as e:
        print(f"Error: Required variable not found in netCDF file: {e}")
        print(f"Available variables: {list(nc.variables.keys())}")
        nc.close()
        raise

    except ValueError as e:
        # This catches the "elevation already exists" error
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
        description="Add elevation variable to a netCDF file based on geopotential (orography)."
    )
    parser.add_argument(
        "filename",
        type=str,
        help="Path to the netCDF file containing geopotential/orography variable"
    )
    parser.add_argument(
        "--geopotential-var",
        type=str,
        default="z",
        help="Name of the geopotential variable (default: z)"
    )
    parser.add_argument(
        "--elevation-var",
        type=str,
        default="elevation",
        help="Name for the output elevation variable (default: elevation)"
    )

    # Parse command-line arguments
    args = parser.parse_args()

    # Process the file
    add_elevation_to_netcdf(
        args.filename,
        geopotential_var=args.geopotential_var,
        elevation_var_name=args.elevation_var
    )
