import argparse
import numpy as np
import netCDF4


def dewpoint_to_specific_humidity(T_d, P):
    """
    Convert dewpoint temperature to specific humidity

    Source:
    -------
    - Stull, R., 2017: "Practical Meteorology: An Algebra-based Survey of Atmospheric
      Science" -version 1.02b.  Univ. of British Columbia.  940 pages.
      isbn 978-0-88865-283-6 .
    - https://www.eoas.ubc.ca/books/Practical_Meteorology/

    Parameters:
    -----------
    T_d : float or array
        Dewpoint temperature [K]
    P : float or array
        Surface pressure [Pa]

    Returns:
    --------
    q : float or array
        Specific humidity [kg/kg]
    """
    # Constants
    epsilon = 0.622  # Ratio of molecular weights

    # Convert dewpoint to vapor pressure using August-Roche-Magnus formula
    # e_s in Pa
    # https://en.wikipedia.org/wiki/Clausius%E2%80%93Clapeyron_relation#Meteorology_and_climatology
    e_s = 610.2 * np.exp(17.625 * (T_d - 273.15) / ((T_d - 273.15) + 243.04))

    # # Tetens formula
    # # Stull2017, Equation (4.2)
    # e_s = 611.3 * np.exp(17.2694 * (T_d - 273.15) / ((T_d - 273.15) + 237.29))

    # Convert vapor pressure to specific humidity
    # Stull2017, Table 4a
    q = epsilon * e_s / (P - (1.0 - epsilon) * e_s)

    return q


def add_specific_humidity_to_netcdf(filename):
    """
    Read surface pressure and dewpoint temperature from a netCDF file,
    calculate specific humidity at 2m, and write it back to the file.

    Parameters:
    -----------
    filename : str
        Path to the netCDF file containing 'sp' and 'd2m' variables

    Returns:
    --------
    None
        Modifies the netCDF file in place by adding 'q2m' variable

    Raises:
    -------
    ValueError
        If 'q2m' variable already exists in the file
    KeyError
        If required variables 'sp' or 'd2m' are not found
    """

    # Open netCDF file in append mode
    print(f"Opening {filename}...")
    nc = netCDF4.Dataset(filename, "a")

    try:
        # Check if q2m already exists - if so, raise error and exit
        if "q2m" in nc.variables:
            nc.close()
            raise ValueError(
                f"Variable 'q2m' already exists in {filename}. "
                "No changes made. Delete the variable first if you want to recalculate."
            )

        # Read the required variables
        print("Reading surface pressure (sp) and dewpoint temperature (d2m)...")
        sp = nc.variables["sp"][:]  # Surface pressure [Pa]
        d2m = nc.variables["d2m"][:]  # Dewpoint temperature at 2m [K]

        print(f"Data shapes - sp: {sp.shape}, d2m: {d2m.shape}")

        # Calculate specific humidity
        print("Calculating specific humidity...")
        q2m = dewpoint_to_specific_humidity(d2m, sp)

        # Create the new variable
        print("Creating new variable 'q2m'...")

        # Get dimension names from d2m (they should be the same)
        dim_names = nc.variables["d2m"].dimensions

        # Create the q2m variable with same dimensions as d2m
        q2m_var = nc.createVariable("q2m", "f4", dim_names, zlib=True, complevel=4)

        # Add attributes
        q2m_var.units = "kg kg-1"
        q2m_var.long_name = "Specific humidity at 2m"
        q2m_var.standard_name = "specific_humidity"
        q2m_var.description = (
            "Calculated from dewpoint temperature (d2m) and surface pressure (sp)"
        )

        # Write the data
        q2m_var[:] = q2m

        print("Variable 'q2m' created and written successfully!")

        # Print some statistics
        print("\nStatistics:")
        print(
            f"  Specific humidity range: {np.min(q2m):.6f} to {np.max(q2m):.6f} kg/kg"
        )
        print(f"  Specific humidity mean: {np.mean(q2m):.6f} kg/kg")
        print(f"  In g/kg: {np.mean(q2m)*1000:.3f} g/kg (mean)")

    except KeyError as e:
        print(f"Error: Required variable not found in netCDF file: {e}")
        print(f"Available variables: {list(nc.variables.keys())}")
        nc.close()
        raise

    except ValueError as e:
        # This catches the "q2m already exists" error
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
        description="Add specific humidity (q2m) to a netCDF file based on dewpoint temperature and surface pressure."
    )
    parser.add_argument(
        "filename",
        type=str,
        help="Path to the ERA5-downloaded netCDF file containing 'sp' and 'd2m' variables"
    )

    # Parse command-line arguments
    args = parser.parse_args()

    # Process the file
    add_specific_humidity_to_netcdf(args.filename)
