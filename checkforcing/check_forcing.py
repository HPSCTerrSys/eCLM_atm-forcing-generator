#!/usr/bin/env python3
"""
Check consistency between eCLM forcing files and namelist configuration.

This script validates that an eCLM directory has consistent atmospheric forcing
data by checking:
    - datm_in namelist file exists and is valid
    - Stream files referenced in datm_in exist
    - Forcing files referenced in stream files exist
    - Forcing files contain required variables

Requirements:
    - f90nml library (pip install f90nml)
    - lxml library (pip install lxml)
    - netCDF4 library (pip install netCDF4)

Usage:
    python check_forcing.py <path_to_datm_in>
    python check_forcing.py /path/to/eclm/run/datm_in
    python check_forcing.py --help
"""

import argparse
# import os
import sys
from pathlib import Path
import f90nml
from lxml import etree
import numpy as np
from netCDF4 import Dataset


def check_forcing_file(forcing_file_path, expected_variables, variable_mappings):
    """Check a NetCDF forcing file for required variables and validate values.

    Args:
        forcing_file_path (Path): Path to the forcing NetCDF file
        expected_variables (list): List of NetCDF variable names expected in file
        variable_mappings (dict): Dict mapping eCLM vars to NetCDF vars

    Returns:
        dict: Dictionary containing check results with keys:
            - 'file_exists': Boolean indicating if file exists
            - 'variables_found': List of found variables
            - 'variables_missing': List of missing variables
            - 'value_checks': Dict of validation results for specific variables
            - 'errors': List of error messages
            - 'warnings': List of warning messages
    """
    result = {
        'file_exists': False,
        'variables_found': [],
        'variables_missing': [],
        'value_checks': {},
        'errors': [],
        'warnings': []
    }

    if not forcing_file_path.exists():
        result['errors'].append(f"File not found: {forcing_file_path}")
        return result

    result['file_exists'] = True

    try:
        with Dataset(str(forcing_file_path), 'r') as nc:
            nc_variables = list(nc.variables.keys())

            # Check which expected variables are present
            for var in expected_variables:
                if var in nc_variables:
                    result['variables_found'].append(var)
                else:
                    result['variables_missing'].append(var)

            # Check for non-negative values in radiation and precipitation variables
            # These eCLM variables should map to non-negative NetCDF variables
            nonneg_eclm_vars = ['swdn', 'lwdn', 'precn']

            for eclm_var in nonneg_eclm_vars:
                if eclm_var in variable_mappings:
                    nc_var = variable_mappings[eclm_var]

                    if nc_var in nc_variables:
                        # Read the variable data
                        var_data = nc.variables[nc_var][:]

                        # Check for negative values
                        min_val = np.min(var_data)
                        max_val = np.max(var_data)
                        has_negative = min_val < 0

                        # Count negative values
                        neg_count = np.sum(var_data < 0)
                        total_count = var_data.size

                        result['value_checks'][nc_var] = {
                            'eclm_var': eclm_var,
                            'min': float(min_val),
                            'max': float(max_val),
                            'has_negative': has_negative,
                            'negative_count': int(neg_count),
                            'total_count': int(total_count)
                        }

                        if has_negative:
                            result['errors'].append(
                                f"{nc_var} (eCLM: {eclm_var}): Found {neg_count} negative "
                                f"values (min={min_val:.4f}), but should be non-negative"
                            )

    except Exception as e:
        result['errors'].append(f"Error reading NetCDF file: {e}")

    return result


def check_stream_forcing_files(root, base_dir):
    """Check all forcing files referenced in a stream file.

    Args:
        root: XML root element from etree.parse()
        base_dir (Path): Base directory for resolving relative paths

    Returns:
        dict: Dictionary with check results for all forcing files
    """
    results = {
        'files_checked': 0,
        'files_valid': 0,
        'files_with_errors': 0,
        'file_results': {},
        'total_errors': 0,
        'total_warnings': 0
    }

    # Extract field info from XML
    field_info = root.find("fieldInfo")
    if field_info is None:
        return results

    # Get forcing file path (can be relative)
    forcing_path_elem = field_info.find("filePath")
    if forcing_path_elem is not None:
        forcing_path = forcing_path_elem.text.strip()
        if not Path(forcing_path).is_absolute():
            forcing_path = base_dir / forcing_path
        else:
            forcing_path = Path(forcing_path)
    else:
        forcing_path = base_dir

    # Get variable mappings (format: "netcdf_var  eclm_var")
    variable_mappings = {}
    var_names = field_info.find("variableNames")
    if var_names is not None:
        var_text = var_names.text.strip()
        for line in var_text.split("\n"):
            line = line.strip()
            if line:
                parts = line.split()
                if len(parts) >= 2:
                    netcdf_var = parts[0]
                    eclm_var = parts[1]
                    variable_mappings[eclm_var] = netcdf_var

    # Get expected NetCDF variable names from mappings
    expected_variables = list(variable_mappings.values())

    # Get forcing file names
    forcing_names = field_info.find("fileNames")
    if forcing_names is None:
        return results

    forcing_files = [line.strip() for line in forcing_names.text.strip().split("\n") if line.strip()]

    # Check each forcing file
    for forcing_file in forcing_files:
        forcing_file_path = forcing_path / forcing_file
        results['files_checked'] += 1

        file_result = check_forcing_file(
            forcing_file_path,
            expected_variables,
            variable_mappings
        )

        results['file_results'][forcing_file] = file_result

        if file_result['errors']:
            results['files_with_errors'] += 1
            results['total_errors'] += len(file_result['errors'])
        else:
            results['files_valid'] += 1

        results['total_warnings'] += len(file_result['warnings'])

    return results


def print_datm_info(nml, base_dir):
    """Print parsed datm_in information.

    Args:
        nml: Namelist object from f90nml.read()
        base_dir (Path): Base directory of datm_in file
    """
    print("\n" + "=" * 70)
    print("DATM Namelist Information")
    print("=" * 70)

    print(f"\nBase directory: {base_dir}")

    # Get domain file if available
    if "shr_strdata_nml" in nml and "domainfile" in nml["shr_strdata_nml"]:
        print(f"Domain file: {nml['shr_strdata_nml']['domainfile']}")

    # Get stream files if available
    if "shr_strdata_nml" in nml and "streams" in nml["shr_strdata_nml"]:
        streams = nml["shr_strdata_nml"]["streams"]
        if isinstance(streams, str):
            streams = [streams]

        # Extract filenames (first part before space)
        stream_files = [s.strip().split()[0] for s in streams if s.strip()]

        print(f"\nStream files found: {len(stream_files)}")
        for i, stream_file in enumerate(stream_files, 1):
            print(f"  {i}. {stream_file}")

    # Print namelist groups
    print("\nNamelist groups:")
    for group_name in nml:
        print(f"  - {group_name}")


def print_stream_info(stream_name, root):
    """Print parsed stream file information.

    Args:
        stream_name (str): Name of the stream file
        root: XML root element from etree.parse()
    """
    print("\n" + "-" * 70)
    print(f"Stream File: {stream_name}")
    print("-" * 70)

    # Get data source
    data_source = root.find("dataSource")
    if data_source is not None:
        print(f"Data source: {data_source.text.strip()}")

    # Get domain info
    domain_info = root.find("domainInfo")
    if domain_info is not None:
        domain_path = domain_info.find("filePath")
        domain_file = domain_info.find("fileNames")

        print("\nDomain:")
        if domain_path is not None:
            print(f"  Path: {domain_path.text.strip()}")
        if domain_file is not None:
            print(f"  File: {domain_file.text.strip()}")

    # Get field info
    field_info = root.find("fieldInfo")
    if field_info is not None:
        forcing_path = field_info.find("filePath")
        forcing_names = field_info.find("fileNames")

        print("\nForcing data:")
        if forcing_path is not None:
            print(f"  Path: {forcing_path.text.strip()}")

        if forcing_names is not None:
            forcing_files = [line.strip() for line in forcing_names.text.strip().split("\n") if line.strip()]
            print(f"  Files: {len(forcing_files)}")

            # Show first few and last few files
            if len(forcing_files) <= 6:
                for f in forcing_files:
                    print(f"    - {f}")
            else:
                for f in forcing_files[:3]:
                    print(f"    - {f}")
                print(f"    ... ({len(forcing_files) - 6} more files)")
                for f in forcing_files[-3:]:
                    print(f"    - {f}")

        # Get variable mappings
        var_names = field_info.find("variableNames")
        if var_names is not None:
            var_text = var_names.text.strip()
            var_mappings = []
            for line in var_text.split("\n"):
                line = line.strip()
                if line:
                    parts = line.split()
                    if len(parts) >= 2:
                        var_mappings.append((parts[0], parts[1]))  # (netcdf_var, eclm_var)

            print(f"\nVariable mappings: {len(var_mappings)}")
            for netcdf_var, eclm_var in var_mappings:
                print(f"  {netcdf_var:12s} (NetCDF) -> {eclm_var} (eCLM)")


def print_forcing_check_results(check_results, verbose=False):
    """Print forcing file check results.

    Args:
        check_results (dict): Output from check_stream_forcing_files()
        verbose (bool): Print detailed results for each file
    """
    print("\nForcing File Checks:")
    print(f"  Files checked: {check_results['files_checked']}")
    print(f"  Files valid: {check_results['files_valid']}")
    print(f"  Files with errors: {check_results['files_with_errors']}")

    if check_results['total_errors'] > 0:
        print(f"  Total errors: {check_results['total_errors']}")

    if check_results['total_warnings'] > 0:
        print(f"  Total warnings: {check_results['total_warnings']}")

    # Show detailed results for files with errors or if verbose
    for filename, file_result in check_results['file_results'].items():
        if file_result['errors'] or verbose:
            print(f"\n  File: {filename}")

            if not file_result['file_exists']:
                print("    ✗ File not found")
                continue

            if file_result['variables_missing']:
                print(f"    ✗ Missing variables: {', '.join(file_result['variables_missing'])}")

            if file_result['variables_found'] and verbose:
                print(f"    ✓ Found variables: {', '.join(file_result['variables_found'])}")

            # Print value check results
            if file_result['value_checks']:
                for nc_var, checks in file_result['value_checks'].items():
                    eclm_var = checks['eclm_var']
                    if checks['has_negative']:
                        print(f"    ✗ {nc_var} (eCLM: {eclm_var}): {checks['negative_count']}/{checks['total_count']} "
                              f"values are negative (min={checks['min']:.4f})")
                    elif verbose:
                        print(f"    ✓ {nc_var} (eCLM: {eclm_var}): All values non-negative "
                              f"(range: [{checks['min']:.4f}, {checks['max']:.4f}])")

            # Print other errors
            for error in file_result['errors']:
                if 'negative values' not in error:  # Already printed above
                    print(f"    ✗ {error}")


def main():
    """Main entry point for the script.

    Returns:
        int: Exit code (0 for success, 1 for error)
    """
    parser = argparse.ArgumentParser(
        description="Parse eCLM datm_in and stream files",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )

    parser.add_argument("datm_in", help="Path to datm_in namelist file")

    parser.add_argument(
        "-v", "--verbose", action="store_true", help="Print verbose output"
    )

    args = parser.parse_args()

    # Parse datm_in using f90nml
    try:
        print(f"Parsing datm_in: {args.datm_in}")
        datm_in_path = Path(args.datm_in).resolve()
        nml = f90nml.read(datm_in_path)
        base_dir = datm_in_path.parent

        print_datm_info(nml, base_dir)
    except Exception as e:
        print(f"\nError parsing datm_in: {e}", file=sys.stderr)
        return 1

    # Extract stream files from namelist
    stream_files = []
    if "shr_strdata_nml" in nml and "streams" in nml["shr_strdata_nml"]:
        streams = nml["shr_strdata_nml"]["streams"]
        if isinstance(streams, str):
            streams = [streams]

        # Extract filenames (first part before space)
        stream_files = [s.strip().split()[0] for s in streams if s.strip()]

    # Parse each stream file and check forcing files
    errors = []
    total_forcing_errors = 0

    for stream_file in stream_files:
        stream_path = base_dir / stream_file

        try:
            # Parse stream XML file using lxml
            tree = etree.parse(str(stream_path))
            root = tree.getroot()

            print_stream_info(stream_file, root)

            # Check forcing files referenced in this stream
            print("\nChecking forcing files...")
            check_results = check_stream_forcing_files(root, base_dir)
            print_forcing_check_results(check_results, verbose=args.verbose)

            total_forcing_errors += check_results['total_errors']

        except FileNotFoundError:
            errors.append(f"Stream file not found: {stream_file}")
            print(f"\n✗ {errors[-1]}", file=sys.stderr)
        except Exception as e:
            errors.append(f"Error parsing {stream_file}: {e}")
            print(f"\n✗ {errors[-1]}", file=sys.stderr)

    # Summary
    print("\n" + "=" * 70)
    print("Summary")
    print("=" * 70)

    if errors:
        print(f"Stream file errors: {len(errors)}")

    if total_forcing_errors > 0:
        print(f"Forcing file errors: {total_forcing_errors}")

    if errors or total_forcing_errors > 0:
        print("\n✗ Validation completed with errors")
        return 1
    else:
        print("✓ All checks passed successfully")
        return 0


if __name__ == "__main__":
    main()
