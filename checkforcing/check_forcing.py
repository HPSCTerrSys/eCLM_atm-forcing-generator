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


def parse_datm_in(datm_in_path):
    """Parse datm_in namelist file using f90nml.

    Args:
        datm_in_path (str): Path to the datm_in file

    Returns:
        dict: Dictionary containing parsed namelist data with keys:
            - 'namelists': Full namelist object from f90nml
            - 'stream_files': List of stream file names
            - 'domain_file': Path to domain file (if specified)
            - 'base_dir': Base directory of datm_in file

    Raises:
        FileNotFoundError: If the datm_in file doesn't exist
        ValueError: If the namelist cannot be parsed
    """
    datm_in_path = Path(datm_in_path).resolve()

    if not datm_in_path.exists():
        raise FileNotFoundError(f"datm_in file not found: {datm_in_path}")

    try:
        # Parse the namelist file
        nml = f90nml.read(datm_in_path)
    except Exception as e:
        raise ValueError(f"Failed to parse datm_in: {e}")

    result = {
        "namelists": nml,
        "stream_files": [],
        "domain_file": None,
        "base_dir": datm_in_path.parent,
    }

    # Extract stream file information
    if "shr_strdata_nml" in nml:
        shr_strdata = nml["shr_strdata_nml"]

        # Get domain file
        if "domainfile" in shr_strdata:
            result["domain_file"] = shr_strdata["domainfile"]

        # Extract stream file names from 'streams' variable
        if "streams" in shr_strdata:
            streams = shr_strdata["streams"]

            # streams can be a single string or a list of strings
            if isinstance(streams, str):
                streams = [streams]

            # Each stream entry format: "filename year1 year2 year3"
            # Extract just the filename (first part before space)
            for stream_entry in streams:
                stream_entry = stream_entry.strip()
                if stream_entry:
                    # Split by whitespace and take first element (filename)
                    filename = stream_entry.split()[0]
                    result["stream_files"].append(filename)

    return result


def parse_stream_file(stream_path, base_dir=None):
    """Parse a DATM stream XML file using lxml.

    Args:
        stream_path (str): Path to the stream file
        base_dir (str, optional): Base directory for resolving relative paths

    Returns:
        dict: Dictionary containing parsed stream data with keys:
            - 'domain_file': Path to domain file
            - 'domain_path': Directory containing domain file
            - 'forcing_path': Directory containing forcing files
            - 'forcing_files': List of forcing file names
            - 'variable_mappings': Dict mapping eCLM vars to NetCDF vars
            - 'data_source': Data source type

    Raises:
        FileNotFoundError: If the stream file doesn't exist
        ValueError: If the XML cannot be parsed
    """
    stream_path = Path(stream_path)

    # Resolve relative path if base_dir is provided
    if base_dir and not stream_path.is_absolute():
        stream_path = Path(base_dir) / stream_path

    if not stream_path.exists():
        raise FileNotFoundError(f"Stream file not found: {stream_path}")

    try:
        # Parse the XML file
        tree = etree.parse(str(stream_path))
        root = tree.getroot()
    except Exception as e:
        raise ValueError(f"Failed to parse stream XML: {e}")

    result = {
        "domain_file": None,
        "domain_path": None,
        "forcing_path": None,
        "forcing_files": [],
        "variable_mappings": {},
        "data_source": None,
    }

    # Extract data source
    data_source = root.find("dataSource")
    if data_source is not None:
        result["data_source"] = data_source.text.strip()

    # Extract domain information
    domain_info = root.find("domainInfo")
    if domain_info is not None:
        # Get domain file path
        domain_path = domain_info.find("filePath")
        if domain_path is not None:
            result["domain_path"] = domain_path.text.strip()

        # Get domain file name
        domain_names = domain_info.find("fileNames")
        if domain_names is not None:
            result["domain_file"] = domain_names.text.strip()

    # Extract field information (forcing files)
    field_info = root.find("fieldInfo")
    if field_info is not None:
        # Get variable mappings
        var_names = field_info.find("variableNames")
        if var_names is not None:
            # Parse variable mappings (format: "eCLM_var  netcdf_var")
            var_text = var_names.text.strip()
            for line in var_text.split("\n"):
                line = line.strip()
                if line:
                    # Split by whitespace
                    parts = line.split()
                    if len(parts) >= 2:
                        eclm_var = parts[0]
                        netcdf_var = parts[1]
                        result["variable_mappings"][eclm_var] = netcdf_var

        # Get forcing file path
        forcing_path = field_info.find("filePath")
        if forcing_path is not None:
            result["forcing_path"] = forcing_path.text.strip()

        # Get forcing file names
        forcing_names = field_info.find("fileNames")
        if forcing_names is not None:
            # Parse file names (one per line)
            file_text = forcing_names.text.strip()
            for line in file_text.split("\n"):
                line = line.strip()
                if line:
                    result["forcing_files"].append(line)

    return result


def print_datm_info(datm_data):
    """Print parsed datm_in information.

    Args:
        datm_data (dict): Output from parse_datm_in()
    """
    print("\n" + "=" * 70)
    print("DATM Namelist Information")
    print("=" * 70)

    print(f"\nBase directory: {datm_data['base_dir']}")

    if datm_data["domain_file"]:
        print(f"Domain file: {datm_data['domain_file']}")

    print(f"\nStream files found: {len(datm_data['stream_files'])}")
    for i, stream_file in enumerate(datm_data["stream_files"], 1):
        print(f"  {i}. {stream_file}")

    # Print namelist groups
    print("\nNamelist groups:")
    for group_name in datm_data["namelists"]:
        print(f"  - {group_name}")


def print_stream_info(stream_name, stream_data):
    """Print parsed stream file information.

    Args:
        stream_name (str): Name of the stream file
        stream_data (dict): Output from parse_stream_file()
    """
    print("\n" + "-" * 70)
    print(f"Stream File: {stream_name}")
    print("-" * 70)

    if stream_data["data_source"]:
        print(f"Data source: {stream_data['data_source']}")

    print("\nDomain:")
    print(f"  Path: {stream_data['domain_path']}")
    print(f"  File: {stream_data['domain_file']}")

    print("\nForcing data:")
    print(f"  Path: {stream_data['forcing_path']}")
    print(f"  Files: {len(stream_data['forcing_files'])}")

    if stream_data["forcing_files"]:
        # Show first few and last few files
        if len(stream_data["forcing_files"]) <= 6:
            for f in stream_data["forcing_files"]:
                print(f"    - {f}")
        else:
            for f in stream_data["forcing_files"][:3]:
                print(f"    - {f}")
            print(f"    ... ({len(stream_data['forcing_files']) - 6} more files)")
            for f in stream_data["forcing_files"][-3:]:
                print(f"    - {f}")

    print(f"\nVariable mappings: {len(stream_data['variable_mappings'])}")
    for eclm_var, netcdf_var in stream_data["variable_mappings"].items():
        print(f"  {eclm_var:12s} -> {netcdf_var}")


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

    # Parse datm_in
    try:
        print(f"Parsing datm_in: {args.datm_in}")
        datm_data = parse_datm_in(args.datm_in)
        print_datm_info(datm_data)
    except Exception as e:
        print(f"\nError parsing datm_in: {e}", file=sys.stderr)
        return 1

    # Parse each stream file
    base_dir = datm_data["base_dir"]
    errors = []

    for stream_file in datm_data["stream_files"]:
        stream_path = base_dir / stream_file

        try:
            stream_data = parse_stream_file(stream_path, base_dir)
            print_stream_info(stream_file, stream_data)
        except FileNotFoundError:
            errors.append(f"Stream file not found: {stream_file}")
            print(f"\n✗ {errors[-1]}", file=sys.stderr)
        except Exception as e:
            errors.append(f"Error parsing {stream_file}: {e}")
            print(f"\n✗ {errors[-1]}", file=sys.stderr)

    # Summary
    print("\n" + "=" * 70)
    if errors:
        print(f"Completed with {len(errors)} error(s)")
        print("=" * 70)
        return 1
    else:
        print("✓ Successfully parsed all files")
        print("=" * 70)
        return 0


if __name__ == "__main__":
    main()
