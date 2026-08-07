#!/usr/bin/env python3
"""
Visualize variables from eCLM atmospheric forcing files.

This script creates time series plots of variables from netCDF forcing files.
Plot labels and metadata are extracted automatically from the netCDF file.

Requirements:
    - xarray
    - matplotlib
    - pandas (for datetime handling)

Usage:
    python plot_forcing.py <file> [<file> ...]
    python plot_forcing.py 2020-09.nc 2020-10.nc --var FSDS --output radiation.png
    python plot_forcing.py --help

Without --var, all time-dependent variables from the first file are plotted.
To list available variables in a file:
    python plot_forcing.py <forcing_file> --list
"""
import argparse
import sys

import matplotlib.pyplot as plt
import xarray as xr


def check_files(file_paths, var_names=None):
    """Sanity checks for multiple files: no time overlap, all variables present.

    Args:
        file_paths (list): List of file paths to check
        var_names (list, optional): Variable names to verify presence in every file
    """
    time_ranges = []
    for path in file_paths:
        try:
            ds = xr.open_dataset(path)
        except FileNotFoundError:
            print(f"Error: File '{path}' not found.")
            sys.exit(1)
        except Exception as e:
            print(f"Error opening file '{path}': {e}")
            sys.exit(1)

        if var_names:
            for vn in var_names:
                if vn not in ds:
                    print(f"Error: Variable '{vn}' not found in '{path}'.")
                    ds.close()
                    sys.exit(1)

        time_ranges.append((ds['time'].values[0], ds['time'].values[-1]))
        ds.close()

    # Check for time overlaps
    sorted_ranges = sorted(time_ranges)
    for i in range(len(sorted_ranges) - 1):
        if sorted_ranges[i][1] >= sorted_ranges[i + 1][0]:
            print("Error: Time overlap detected between input files.")
            sys.exit(1)


def list_variables(ds):
    """List all data variables in the dataset with their descriptions.

    Args:
        ds (xr.Dataset): The xarray dataset

    Returns:
        list: List of (variable_name, long_name, units) tuples
    """
    variables = []
    for var_name in ds.data_vars:
        var = ds[var_name]
        # Skip coordinate variables
        if 'time' not in var.dims:
            continue
        long_name = var.attrs.get('long_name', var_name)
        units = var.attrs.get('units', '')
        variables.append((var_name, long_name, units))
    return variables


def get_location_info(ds):
    """Extract location information from dataset coordinates.

    Args:
        ds (xr.Dataset): The xarray dataset

    Returns:
        tuple: (latitude, longitude) or (None, None) if not found
    """
    lat, lon = None, None

    # Try xc/yc coordinates
    if 'xc' in ds and 'yc' in ds:
        lon = float(ds['xc'].values.flatten()[0])
        lat = float(ds['yc'].values.flatten()[0])

    return lat, lon


def plot_variable(ds, var_name, ax=None, output_path=None, show=True):
    """Create a time series plot for a variable.

    Args:
        ds (xr.Dataset): The xarray dataset
        var_name (str): Name of the variable to plot
        ax (matplotlib.axes.Axes, optional): Axes to plot into. If None, a new
            figure is created and output_path/show are handled here.
        output_path (str, optional): Path to save the figure (standalone only)
        show (bool): Whether to display the plot interactively (standalone only)

    Returns:
        matplotlib.figure.Figure: The created figure
    """
    if var_name not in ds:
        raise ValueError(f"Variable '{var_name}' not found in dataset. "
                         f"Available variables: {list(ds.data_vars)}")

    var = ds[var_name]

    # Check if variable has time dimension
    if 'time' not in var.dims:
        raise ValueError(f"Variable '{var_name}' does not have a time dimension. "
                         "Only time-dependent variables can be plotted.")

    # Extract metadata
    long_name = var.attrs.get('long_name', var_name)
    units = var.attrs.get('units', '')
    institution = var.attrs.get('institution', '')

    # Get location info
    lat, lon = get_location_info(ds)

    # Squeeze out singleton dimensions (nj, ni)
    data = var.squeeze()

    # Create figure if no axes provided
    standalone = ax is None
    if standalone:
        fig, ax = plt.subplots(figsize=(12, 5))
    else:
        fig = ax.get_figure()

    # Plot the data
    data.plot(ax=ax, linewidth=0.8)

    # Set labels
    ylabel = f"{var_name}"
    if units:
        ylabel += f" [{units}]"
    ax.set_ylabel(ylabel)
    ax.set_xlabel("Time")

    # Create title
    title = f"{var_name}: {long_name}"
    if lat is not None and lon is not None:
        title += f"\nLocation: {lat:.3f}N, {lon:.3f}E"
    if institution:
        title += f" (Source: {institution})"
    ax.set_title(title)

    ax.grid(True, alpha=0.3)

    if standalone:
        fig.tight_layout()
        if output_path:
            fig.savefig(output_path, dpi=150, bbox_inches='tight')
            print(f"Figure saved to: {output_path}")
        if show:
            plt.show()

    return fig


def main():
    """Main function to parse arguments and create plots."""
    parser = argparse.ArgumentParser(
        description="Visualize variables from eCLM atmospheric forcing files.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
    python plot_forcing.py 2020-09.nc
    python plot_forcing.py 2020-09.nc --var TBOT FSDS
    python plot_forcing.py 2020-09.nc 2020-10.nc --var FSDS --output radiation.png
    python plot_forcing.py 2020-09.nc --list
        """
    )
    parser.add_argument(
        "forcing_files",
        type=str,
        nargs='+',
        help="Path(s) to the forcing netCDF file(s)"
    )
    parser.add_argument(
        "--var", "-v",
        type=str,
        nargs='+',
        default=None,
        help="Variable(s) to plot (e.g., TBOT FSDS). Plots all variables if omitted."
    )
    parser.add_argument(
        "--list", "-l",
        action="store_true",
        help="List all available time-dependent variables in the file"
    )
    parser.add_argument(
        "--output", "-o",
        type=str,
        default=None,
        help="Output file path for the figure (only used with --var)"
    )
    parser.add_argument(
        "--no-show",
        action="store_true",
        help="Do not display the plot interactively (useful for batch processing)"
    )

    args = parser.parse_args()

    # Open first file to determine variable list
    try:
        ds0 = xr.open_dataset(args.forcing_files[0])
    except FileNotFoundError:
        print(f"Error: File '{args.forcing_files[0]}' not found.")
        sys.exit(1)
    except Exception as e:
        print(f"Error opening file: {e}")
        sys.exit(1)

    # List variables if requested
    if args.list:
        variables = list_variables(ds0)
        print(f"\nTime-dependent variables in '{args.forcing_files[0]}':")
        print("-" * 70)
        print(f"{'Variable':<12} {'Long Name':<40} {'Units':<15}")
        print("-" * 70)
        for var_name, long_name, units in variables:
            if len(long_name) > 38:
                long_name = long_name[:35] + "..."
            print(f"{var_name:<12} {long_name:<40} {units:<15}")
        print("-" * 70)
        ds0.close()
        return

    # Determine variables to plot
    if args.var:
        var_names = args.var
    else:
        var_names = [vn for vn, _, _ in list_variables(ds0)]
        if not var_names:
            print("Error: No time-dependent variables found in the first file.")
            ds0.close()
            sys.exit(1)
    ds0.close()

    # Sanity checks and open dataset(s)
    if len(args.forcing_files) > 1:
        check_files(args.forcing_files, var_names=var_names)
        try:
            ds = xr.open_mfdataset(args.forcing_files, combine='by_coords')
        except Exception as e:
            print(f"Error opening files: {e}")
            sys.exit(1)
    else:
        try:
            ds = xr.open_dataset(args.forcing_files[0])
        except FileNotFoundError:
            print(f"Error: File '{args.forcing_files[0]}' not found.")
            sys.exit(1)
        except Exception as e:
            print(f"Error opening file: {e}")
            sys.exit(1)

    # Create the plot(s)
    try:
        if len(var_names) == 1:
            plot_variable(
                ds, var_names[0],
                output_path=args.output,
                show=not args.no_show
            )
        else:
            fig, axes = plt.subplots(
                len(var_names), 1,
                figsize=(12, 4 * len(var_names)),
                sharex=True
            )
            for var_name, ax in zip(var_names, axes):
                plot_variable(ds, var_name, ax=ax)
            fig.tight_layout()
            if args.output:
                fig.savefig(args.output, dpi=150, bbox_inches='tight')
                print(f"Figure saved to: {args.output}")
            if not args.no_show:
                plt.show()
    except ValueError as e:
        print(f"Error: {e}")
        sys.exit(1)
    finally:
        ds.close()


if __name__ == "__main__":
    main()
