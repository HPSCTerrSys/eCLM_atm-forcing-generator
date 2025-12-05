import argparse
import numpy as np
import netCDF4
# import xarray as xr
from typing import Union


def temperature_extrapolation_adiabatic(
    T_2m: Union[float, np.ndarray],
    z1: float = 2.0,
    z2: float = 10.0,
    lapse_rate: float = -0.0065,
) -> Union[float, np.ndarray]:
    """
    Temperature extrapolation using atmospheric lapse rate.

    Formula:
    --------
    T(z2) = T(z1) + Γ * (z2 - z1)

    where Γ is the atmospheric lapse rate (K/m)

    Parameters:
    -----------
    T_2m : float or array
        Temperature at 2m height (K or °C)
    z1 : float
        Initial height (m), default 2.0
    z2 : float
        Target height (m), default 10.0
    lapse_rate : float
        Temperature lapse rate (K/m), default
        -0.0065 K/m (free atmosphere)

    Returns:
    --------
    T_10m : float or array
        Temperature at 10m height (same units as input)

    References:
    -----------
    Lapse rate from standard atmosphere (ISO 2533:1975).
    """
    dz = z2 - z1
    T_10m = T_2m + lapse_rate * dz
    return T_10m


def humidity_extrapolation_constant_mixing_ratio(
    q_2m: Union[float, np.ndarray],
) -> Union[float, np.ndarray]:
    """
    Extrapolate specific humidity assuming constant mixing ratio.

    The mixing ratio is assumed constant with height.
    In this simple case: q(2m) = q(10m).

    Formula:
    --------
    1. Calculate mixing ratio at 2m:
       r = q / (1 - q)  [kg/kg]

    2. Assume mixing ration at 10m:
       r(10m) = r(2m)

    3. Convert back to specific humidity:
       q(10m) = r / (1 + r)

    Parameters:
    -----------
    q_2m : float or array
        Specific humidity at 2m (kg/kg)

    Returns:
    --------
    q_10m : float or array
        Specific humidity at 10m (kg/kg)

    """
    # Convert to mixing ratio
    mixing_ratio_2m = q_2m / (1.0 - q_2m)

    # Assume constant mixing ratio (well-mixed assumption)
    mixing_ratio_10m = mixing_ratio_2m

    # Convert back to specific humidity
    q_10m = mixing_ratio_10m / (1.0 + mixing_ratio_10m)

    return q_10m


def convert_2m_to_10m_in_netcdf(filename):
    """
    Read 2m temperature and specific humidity from a netCDF file,
    extrapolate to 10m height, and write back to the file.

    Parameters:
    -----------
    filename : str
        Path to the netCDF file containing 't2m' and 'q2m' variables

    Returns:
    --------
    None
        Modifies the netCDF file in place by adding 't10m' and 'q10m' variables

    Raises:
    -------
    ValueError
        If 't10m' or 'q10m' variables already exist in the file
    KeyError
        If required variables 't2m' or 'q2m' are not found
    """

    # Open netCDF file in append mode
    print(f"Opening {filename}...")
    nc = netCDF4.Dataset(filename, "a")

    try:
        # Check if t10m or q10m already exist - if so, raise error and exit
        if "t10m" in nc.variables or "q10m" in nc.variables:
            nc.close()
            raise ValueError(
                f"Variable 't10m' and/or 'q10m' already exist in {filename}. "
                "No changes made. Delete the variable(s) first if you want to recalculate."
            )

        # Read the required variables
        print("Reading temperature (t2m) and specific humidity (q2m)...")
        t2m = nc.variables["t2m"][:]  # Temperature at 2m [K]
        q2m = nc.variables["q2m"][:]  # Specific humidity at 2m [kg/kg]

        print(f"Data shapes - t2m: {t2m.shape}, q2m: {q2m.shape}")

        # Extrapolate to 10m
        print("Extrapolating temperature to 10m using adiabatic lapse rate...")
        t10m = temperature_extrapolation_adiabatic(t2m)

        print(
            "Extrapolating specific humidity to 10m assuming constant mixing ratio..."
        )
        q10m = humidity_extrapolation_constant_mixing_ratio(q2m)

        # Create the new variables
        print("Creating new variables 't10m' and 'q10m'...")

        # Get dimension names from t2m
        dim_names = nc.variables["t2m"].dimensions

        # Create the t10m variable
        t10m_var = nc.createVariable("t10m", "f4", dim_names, zlib=True, complevel=4)

        # Add attributes for t10m
        t10m_var.units = "K"
        t10m_var.long_name = "Temperature at 10m"
        t10m_var.standard_name = "air_temperature"
        t10m_var.description = (
            "Extrapolated from 2m temperature (t2m) using adiabatic lapse rate"
        )

        # Write the temperature data
        t10m_var[:] = t10m

        # Create the q10m variable
        q10m_var = nc.createVariable("q10m", "f4", dim_names, zlib=True, complevel=4)

        # Add attributes for q10m
        q10m_var.units = "kg kg-1"
        q10m_var.long_name = "Specific humidity at 10m"
        q10m_var.standard_name = "specific_humidity"
        q10m_var.description = "Extrapolated from 2m specific humidity (q2m) assuming constant mixing ratio"

        # Write the humidity data
        q10m_var[:] = q10m

        print("Variables 't10m' and 'q10m' created and written successfully!")

        # Print some statistics
        print("\nStatistics:")
        print(f"  Temperature 2m range: {np.min(t2m):.2f} to {np.max(t2m):.2f} K")
        print(f"  Temperature 10m range: {np.min(t10m):.2f} to {np.max(t10m):.2f} K")
        print(f"  Temperature difference (10m-2m) mean: {np.mean(t10m - t2m):.4f} K")
        print(f"  Specific humidity 2m mean: {np.mean(q2m):.6f} kg/kg")
        print(f"  Specific humidity 10m mean: {np.mean(q10m):.6f} kg/kg")

    except KeyError as e:
        print(f"Error: Required variable not found in netCDF file: {e}")
        print(f"Available variables: {list(nc.variables.keys())}")
        nc.close()
        raise

    except ValueError as e:
        # This catches the "t10m/q10m already exists" error
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
        description="Convert 2m temperature and humidity to 10m height in a netCDF file."
    )
    parser.add_argument(
        "filename",
        type=str,
        help="Path to the ERA5-downloaded netCDF file containing 't2m' and 'q2m' variables",
    )

    # Parse command-line arguments
    args = parser.parse_args()

    # Process the file
    convert_2m_to_10m_in_netcdf(args.filename)
