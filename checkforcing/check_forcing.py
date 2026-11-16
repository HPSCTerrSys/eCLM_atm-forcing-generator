#!/usr/bin/env python3
"""
Check consistency and validity of eCLM atmospheric forcing files.

Two input modes:
  Namelist mode:  check_forcing.py <datm_in>
                  Full validation chain: namelist -> stream XML -> NetCDF files.
  Directory mode: check_forcing.py --dir <forcing_dir>
                  Directly scan NetCDF files; checks variable presence and
                  non-negative values (no namelist/stream validation).

Ensemble mode (--ensemble-size N):
  Namelist mode:  processes datm_in_0001 .. datm_in_NNNN as siblings of datm_in
  Directory mode: processes real_00001 .. real_NNNNN under --dir

Requirements:
    - f90nml library (pip install f90nml)   [namelist mode only]
    - lxml library (pip install lxml)       [namelist mode only]
    - netCDF4 library (pip install netCDF4)
    - numpy library (pip install numpy)

Usage:
    python check_forcing.py datm_in
    python check_forcing.py datm_in --ensemble-size 50
    python check_forcing.py --dir /path/to/forcings
    python check_forcing.py --dir /path/to/run --ensemble-size 50 --members 1 2 3
    python check_forcing.py --help
"""

import argparse
import sys
from pathlib import Path
import numpy as np
from netCDF4 import Dataset


# ---------------------------------------------------------------------------
# Constants for directory mode
# ---------------------------------------------------------------------------

DEFAULT_DIR_VARIABLES = ["FSDS", "FLDS", "PRECTmms", "TBOT", "WIND", "PSRF", "QBOT"]
# NetCDF variable name -> eCLM internal name (triggers non-negative check in check_forcing_file)
_NONNEG_NC_TO_ECLM = {"FSDS": "swdn", "FLDS": "lwdn", "PRECTmms": "precn"}


def _build_dir_variable_mappings(variables):
    """Build variable_mappings dict for directory mode.

    Non-negative variables get their eCLM alias as key so check_forcing_file()
    triggers the value check. All others use an identity mapping (presence only).
    """
    mappings = {}
    for v in variables:
        if v in _NONNEG_NC_TO_ECLM:
            mappings[_NONNEG_NC_TO_ECLM[v]] = v  # e.g. "swdn" -> "FSDS"
        else:
            mappings[v] = v
    return mappings


# ---------------------------------------------------------------------------
# Per-file check and fix
# ---------------------------------------------------------------------------

def check_forcing_file(forcing_file_path, variable_mappings):
    """Check a NetCDF forcing file for required variables and validate values.

    Args:
        forcing_file_path (Path): Path to the forcing NetCDF file
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
        "file_exists": False,
        "variables_found": [],
        "variables_missing": [],
        "value_checks": {},
        "errors": [],
        "warnings": [],
    }

    if not forcing_file_path.exists():
        result["errors"].append(f"File not found: {forcing_file_path}")
        return result

    result["file_exists"] = True

    try:
        with Dataset(str(forcing_file_path), "r") as nc:
            nc_variables = list(nc.variables.keys())

            # Check which expected netCDF variables are present in forcing file
            for var in variable_mappings.values():
                if var in nc_variables:
                    result["variables_found"].append(var)
                else:
                    result["variables_missing"].append(var)

            # Check for non-negative values in radiation and precipitation variables
            # These eCLM variables should map to non-negative NetCDF variables
            nonneg_eclm_vars = ["swdn", "lwdn", "precn"]

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

                        result["value_checks"][nc_var] = {
                            "eclm_var": eclm_var,
                            "min": float(min_val),
                            "max": float(max_val),
                            "has_negative": has_negative,
                            "negative_count": int(neg_count),
                            "total_count": int(total_count),
                        }

                        if has_negative:
                            result["errors"].append(
                                f"{nc_var} (eCLM: {eclm_var}): Found {neg_count} negative "
                                f"values (min={min_val:.4f}), but should be non-negative"
                            )

    except Exception as e:
        result["errors"].append(f"Error reading forcing file: {e}")

    return result


def fix_negative_values_in_file(forcing_file_path, variable_mappings):
    """Fix negative values in a NetCDF forcing file by setting them to zero.

    Args:
        forcing_file_path (Path): Path to the forcing NetCDF file
        variable_mappings (dict): Dict mapping eCLM vars to NetCDF vars

    Returns:
        dict: Dictionary with fix results:
            - 'fixed': Boolean indicating if any values were fixed
            - 'variables_fixed': Dict of variable names and count of values fixed
            - 'error': Error message if fixing failed
    """
    result = {
        'fixed': False,
        'variables_fixed': {},
        'error': None
    }

    if not forcing_file_path.exists():
        result['error'] = f"File not found: {forcing_file_path}"
        return result

    try:
        # Open in write mode
        with Dataset(str(forcing_file_path), 'r+') as nc:
            nc_variables = list(nc.variables.keys())

            # Variables that should be non-negative
            nonneg_eclm_vars = ['swdn', 'lwdn', 'precn']
            fixed_vars = []

            for eclm_var in nonneg_eclm_vars:
                if eclm_var in variable_mappings:
                    nc_var = variable_mappings[eclm_var]

                    if nc_var in nc_variables:
                        # Read the variable data
                        var_data = nc.variables[nc_var][:]

                        # Check for negative values
                        neg_mask = var_data < 0
                        neg_count = np.sum(neg_mask)

                        if neg_count > 0:
                            # Set negative values to zero
                            var_data[neg_mask] = 0
                            nc.variables[nc_var][:] = var_data

                            result['variables_fixed'][nc_var] = int(neg_count)
                            fixed_vars.append(f"{nc_var}({neg_count})")
                            result['fixed'] = True

            # Update history attribute if we fixed anything
            if result['fixed']:
                import datetime
                timestamp = datetime.datetime.now().strftime("%a %b %d %H:%M:%S %Y")
                history_entry = (
                    f"{timestamp}: check_forcing.py --fix-negative: "
                    f"Set negative values to zero in {', '.join(fixed_vars)}\n"
                )

                # Prepend to existing history or create new one
                if hasattr(nc, 'history'):
                    nc.history = history_entry + nc.history
                else:
                    nc.history = history_entry.strip()

    except Exception as e:
        result['error'] = f"Error fixing file: {e}"

    return result


# ---------------------------------------------------------------------------
# Stream-based check
# ---------------------------------------------------------------------------

def check_stream_forcing_files(root, base_dir, fix_negative=False):
    """Check all forcing files referenced in a stream file.

    Args:
        root: XML root element from etree.parse()
        base_dir (Path): Base directory for resolving relative paths
        fix_negative (bool): Fix negative values by setting them to zero

    Returns:
        dict: Dictionary with check results for all forcing files
    """
    results = {
        "files_checked": 0,
        "files_valid": 0,
        "files_with_errors": 0,
        "files_fixed": 0,
        "file_results": {},
        "total_errors": 0,
        "total_warnings": 0,
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

    # Get forcing file names
    forcing_names = field_info.find("fileNames")
    if forcing_names is None:
        return results

    forcing_files = [
        line.strip() for line in forcing_names.text.strip().split("\n") if line.strip()
    ]

    # Check each forcing file
    for forcing_file in forcing_files:
        forcing_file_path = forcing_path / forcing_file
        results["files_checked"] += 1

        file_result = check_forcing_file(forcing_file_path, variable_mappings)

        # Check if file has negative value errors
        has_negative_errors = False
        if file_result["value_checks"]:
            for checks in file_result["value_checks"].values():
                if checks["has_negative"]:
                    has_negative_errors = True
                    break

        # Fix negative values if requested and errors exist
        if fix_negative and has_negative_errors:
            fix_result = fix_negative_values_in_file(forcing_file_path, variable_mappings)
            file_result["fix_result"] = fix_result

            if fix_result["fixed"]:
                results["files_fixed"] += 1
                # Re-check the file after fixing
                file_result = check_forcing_file(forcing_file_path, variable_mappings)
                file_result["was_fixed"] = True
                file_result["fix_result"] = fix_result

        results["file_results"][forcing_file] = file_result

        if file_result["errors"]:
            results["files_with_errors"] += 1
            results["total_errors"] += len(file_result["errors"])
        else:
            results["files_valid"] += 1

        results["total_warnings"] += len(file_result["warnings"])

    return results


# ---------------------------------------------------------------------------
# Directory-mode check
# ---------------------------------------------------------------------------

def check_member_directory(dir_path, variables=None, fix_negative=False, verbose=False):
    """Check all NetCDF files in dir_path directly (no namelist/stream parsing).

    Uses the same result structure as check_stream_forcing_files so that
    print_forcing_check_results can be reused unchanged.

    Args:
        dir_path (Path): Directory containing *.nc forcing files
        variables (list[str]): Variable names to check. Defaults to DEFAULT_DIR_VARIABLES.
        fix_negative (bool): Fix negative values by setting them to zero
        verbose (bool): Print verbose output

    Returns:
        dict: Same structure as check_stream_forcing_files output
    """
    if variables is None:
        variables = DEFAULT_DIR_VARIABLES

    variable_mappings = _build_dir_variable_mappings(variables)

    results = {
        "files_checked": 0,
        "files_valid": 0,
        "files_with_errors": 0,
        "files_fixed": 0,
        "file_results": {},
        "total_errors": 0,
        "total_warnings": 0,
    }

    nc_files = sorted(dir_path.glob("*.nc"))
    if not nc_files:
        print(f"  (no .nc files found in {dir_path})")
        return results

    if verbose:
        print(f"\nChecking {len(nc_files)} NetCDF file(s) in {dir_path}")
        print(f"Variables: {', '.join(variables)}")

    for nc_file in nc_files:
        results["files_checked"] += 1

        file_result = check_forcing_file(nc_file, variable_mappings)

        has_negative_errors = any(
            checks["has_negative"]
            for checks in file_result["value_checks"].values()
        )

        if fix_negative and has_negative_errors:
            fix_result = fix_negative_values_in_file(nc_file, variable_mappings)
            file_result["fix_result"] = fix_result

            if fix_result["fixed"]:
                results["files_fixed"] += 1
                file_result = check_forcing_file(nc_file, variable_mappings)
                file_result["was_fixed"] = True
                file_result["fix_result"] = fix_result

        results["file_results"][nc_file.name] = file_result

        if file_result["errors"]:
            results["files_with_errors"] += 1
            results["total_errors"] += len(file_result["errors"])
        else:
            results["files_valid"] += 1

        results["total_warnings"] += len(file_result["warnings"])

    return results


# ---------------------------------------------------------------------------
# Namelist-mode check (single member)
# ---------------------------------------------------------------------------

def check_member_namelist(datm_in_path, fix_negative=False, verbose=False):
    """Run the full namelist-mode check for a single datm_in file.

    Prints per-stream results directly. Returns a summary dict for use by
    main() (single mode) and ensemble aggregation.

    Args:
        datm_in_path (Path): Path to the datm_in namelist file
        fix_negative (bool): Fix negative values by setting them to zero
        verbose (bool): Print verbose output

    Returns:
        dict with keys:
            - 'ok': bool
            - 'total_errors': int  (forcing file errors)
            - 'stream_errors': int (stream file parse/not-found errors)
            - 'fsds_min': float | None
    """
    import f90nml
    from lxml import etree

    summary = {"ok": True, "total_errors": 0, "stream_errors": 0, "fsds_min": None}

    # Parse datm_in
    try:
        if verbose:
            print(f"Parsing datm_in: {datm_in_path}")
        nml = f90nml.read(datm_in_path)
        base_dir = datm_in_path.parent
        print_datm_info(nml, base_dir, verbose=verbose)
    except Exception as e:
        print(f"\nError parsing datm_in: {e}", file=sys.stderr)
        summary["ok"] = False
        summary["stream_errors"] += 1
        return summary

    # Extract stream files from namelist
    stream_files = []
    if "shr_strdata_nml" in nml and "streams" in nml["shr_strdata_nml"]:
        streams = nml["shr_strdata_nml"]["streams"]
        if isinstance(streams, str):
            streams = [streams]
        stream_files = [s.strip().split()[0] for s in streams if s.strip()]

    all_stream_results = []

    for stream_file in stream_files:
        stream_path = base_dir / stream_file

        try:
            tree = etree.parse(str(stream_path))
            root = tree.getroot()

            print_stream_info(stream_file, root, verbose=verbose)

            if verbose:
                print("\nChecking forcing files...")
            check_results = check_stream_forcing_files(
                root, base_dir, fix_negative=fix_negative
            )
            print_forcing_check_results(check_results, verbose=verbose)

            summary["total_errors"] += check_results["total_errors"]
            all_stream_results.append(check_results)

        except FileNotFoundError:
            print(f"✗ Stream file not found: {stream_file}", file=sys.stderr)
            summary["stream_errors"] += 1
        except Exception as e:
            print(f"✗ Error parsing {stream_file}: {e}", file=sys.stderr)
            summary["stream_errors"] += 1

    # Extract FSDS min across all stream results for ensemble summary
    fsds_mins = [
        file_result["value_checks"]["FSDS"]["min"]
        for sr in all_stream_results
        for file_result in sr["file_results"].values()
        if "FSDS" in file_result.get("value_checks", {})
    ]
    if fsds_mins:
        summary["fsds_min"] = min(fsds_mins)

    if summary["total_errors"] > 0 or summary["stream_errors"] > 0:
        summary["ok"] = False

    return summary


# ---------------------------------------------------------------------------
# Print helpers
# ---------------------------------------------------------------------------

def print_datm_info(nml, base_dir, verbose=False):
    """Print parsed datm_in information.

    Args:
        nml: Namelist object from f90nml.read()
        base_dir (Path): Base directory of datm_in file
        verbose (bool): Print detailed output
    """
    if not verbose:
        return

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


def print_stream_info(stream_name, root, verbose=False):
    """Print parsed stream file information.

    Args:
        stream_name (str): Name of the stream file
        root: XML root element from etree.parse()
        verbose (bool): Print detailed output
    """
    if not verbose:
        return

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
            forcing_files = [
                line.strip()
                for line in forcing_names.text.strip().split("\n")
                if line.strip()
            ]
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
                        var_mappings.append(
                            (parts[0], parts[1])
                        )  # (netcdf_var, eclm_var)

            print(f"\nVariable mappings: {len(var_mappings)}")
            for netcdf_var, eclm_var in var_mappings:
                print(f"  {netcdf_var:12s} (NetCDF) -> {eclm_var} (eCLM)")


def print_forcing_check_results(check_results, verbose=False):
    """Print forcing file check results.

    Args:
        check_results (dict): Output from check_stream_forcing_files() or
                              check_member_directory()
        verbose (bool): Print detailed results for each file
    """
    # In non-verbose mode, only show errors and fixes
    if not verbose:
        # Show files that were fixed
        for filename, file_result in check_results["file_results"].items():
            if file_result.get("was_fixed"):
                fix_info = file_result.get("fix_result", {})
                if fix_info.get("fixed"):
                    vars_fixed = ", ".join([
                        f"{var}({count})"
                        for var, count in fix_info["variables_fixed"].items()
                    ])
                    print(f"\n  ✓ Fixed {filename}: {vars_fixed}")

        # Show files with errors (that weren't fixed or still have errors)
        for filename, file_result in check_results["file_results"].items():
            if file_result["errors"] and not file_result.get("was_fixed"):
                print(f"\n  ✗ {filename}:")

                if not file_result["file_exists"]:
                    print("      File not found")
                    continue

                if file_result["variables_missing"]:
                    print(
                        f"      Missing variables: {', '.join(file_result['variables_missing'])}"
                    )

                # Print value check results (errors only)
                if file_result["value_checks"]:
                    for nc_var, checks in file_result["value_checks"].items():
                        if checks["has_negative"]:
                            eclm_var = checks["eclm_var"]
                            print(
                                f"      {nc_var} (eCLM: {eclm_var}): {checks['negative_count']}/{checks['total_count']} "
                                f"values are negative (min={checks['min']:.4f})"
                            )

                # Print other errors
                for error in file_result["errors"]:
                    if "negative values" not in error:  # Already printed above
                        print(f"      {error}")
        return

    # Verbose mode: show everything
    print("\nForcing File Checks:")
    print(f"  Files checked: {check_results['files_checked']}")
    print(f"  Files valid: {check_results['files_valid']}")
    print(f"  Files with errors: {check_results['files_with_errors']}")

    if check_results.get("files_fixed", 0) > 0:
        print(f"  Files fixed: {check_results['files_fixed']}")

    if check_results["total_errors"] > 0:
        print(f"  Total errors: {check_results['total_errors']}")

    if check_results["total_warnings"] > 0:
        print(f"  Total warnings: {check_results['total_warnings']}")

    # Show detailed results for files with errors or if verbose
    for filename, file_result in check_results["file_results"].items():
        if file_result["errors"] or verbose:
            print(f"\n  File: {filename}")

            # Show fix info if file was fixed
            if file_result.get("was_fixed"):
                fix_info = file_result.get("fix_result", {})
                if fix_info.get("fixed"):
                    vars_fixed = ", ".join([
                        f"{var}({count})"
                        for var, count in fix_info["variables_fixed"].items()
                    ])
                    print(f"    ✓ Fixed negative values: {vars_fixed}")

            if not file_result["file_exists"]:
                print("    ✗ File not found")
                continue

            if file_result["variables_missing"]:
                print(
                    f"    ✗ Missing variables: {', '.join(file_result['variables_missing'])}"
                )

            if file_result["variables_found"] and verbose:
                print(
                    f"    ✓ Found variables: {', '.join(file_result['variables_found'])}"
                )

            # Print value check results
            if file_result["value_checks"]:
                for nc_var, checks in file_result["value_checks"].items():
                    eclm_var = checks["eclm_var"]
                    if checks["has_negative"]:
                        print(
                            f"    ✗ {nc_var} (eCLM: {eclm_var}): {checks['negative_count']}/{checks['total_count']} "
                            f"values are negative (min={checks['min']:.4f})"
                        )
                    elif verbose:
                        print(
                            f"    ✓ {nc_var} (eCLM: {eclm_var}): All values non-negative "
                            f"(range: [{checks['min']:.4f}, {checks['max']:.4f}])"
                        )

            # Print other errors
            for error in file_result["errors"]:
                if "negative values" not in error:  # Already printed above
                    print(f"    ✗ {error}")


def print_ensemble_summary(all_results):
    """Print a cross-member summary table.

    Args:
        all_results (list[dict]): Each dict has keys: member_id, ok,
                                  total_errors, stream_errors, fsds_min
    """
    print("\n" + "=" * 70)
    print("Ensemble Summary")
    print("=" * 70)
    print(f"  {'Member':<20}  {'Status':<6}  {'Errors':>6}  {'FSDS min':>10}")
    print(f"  {'-'*20}  {'-'*6}  {'-'*6}  {'-'*10}")

    n_ok = 0
    n_fail = 0
    for r in all_results:
        status = "OK" if r["ok"] else "FAIL"
        errors = r["total_errors"] + r.get("stream_errors", 0)
        fsds_str = f"{r['fsds_min']:.4f}" if r["fsds_min"] is not None else "N/A"
        print(f"  {r['member_id']:<20}  {status:<6}  {errors:>6}  {fsds_str:>10}")
        if r["ok"]:
            n_ok += 1
        else:
            n_fail += 1

    print()
    print(f"  Total: {len(all_results)} member(s)  |  {n_ok} OK  |  {n_fail} FAIL")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    """Main entry point for the script.

    Returns:
        int: Exit code (0 for success, 1 for error)
    """
    parser = argparse.ArgumentParser(
        description="Check eCLM atmospheric forcing files",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )

    parser.add_argument(
        "datm_in", nargs="?", default=None,
        help="Path to datm_in namelist file (namelist mode)",
    )
    parser.add_argument(
        "--dir", metavar="FORCING_DIR",
        help="Forcing directory to scan directly (directory mode; "
             "mutually exclusive with datm_in)",
    )
    parser.add_argument(
        "--ensemble-size", type=int, default=None, metavar="N",
        help="Number of ensemble members. Namelist mode: looks for datm_in_NNNN "
             "siblings of datm_in. Directory mode: looks for real_NNNNN subdirs "
             "under --dir.",
    )
    parser.add_argument(
        "--members", nargs="+", type=int, default=None, metavar="I",
        help="Member indices to check (1-based). Default: all 1..ensemble-size.",
    )
    parser.add_argument(
        "--variables", nargs="+", default=None, metavar="VAR",
        help=f"Variables to check in directory mode "
             f"(default: {' '.join(DEFAULT_DIR_VARIABLES)})",
    )
    parser.add_argument(
        "-v", "--verbose", action="store_true", help="Print verbose output",
    )
    parser.add_argument(
        "--fix-negative", action="store_true",
        help="Fix negative values by setting them to zero and update file history",
    )

    args = parser.parse_args()

    # Validate input mode
    if args.datm_in and args.dir:
        parser.error("Specify either datm_in or --dir, not both.")
    if not args.datm_in and not args.dir:
        parser.error("Specify either datm_in (namelist mode) or --dir (directory mode).")

    dir_mode = args.dir is not None
    variables = args.variables or DEFAULT_DIR_VARIABLES

    # Determine member list
    if args.ensemble_size:
        members = args.members if args.members else list(range(1, args.ensemble_size + 1))
    else:
        members = None

    # -----------------------------------------------------------------------
    # Single mode
    # -----------------------------------------------------------------------
    if members is None:
        if dir_mode:
            dir_path = Path(args.dir).resolve()
            if not dir_path.exists():
                sys.exit(f"ERROR: directory '{dir_path}' not found.")
            print(f"Checking directory: {dir_path}")
            check_results = check_member_directory(
                dir_path, variables=variables,
                fix_negative=args.fix_negative, verbose=args.verbose,
            )
            print_forcing_check_results(check_results, verbose=args.verbose)

            print("\n" + "=" * 70)
            print("Summary")
            print("=" * 70)
            print(f"Files checked: {check_results['files_checked']}")
            if check_results["total_errors"] > 0:
                print(f"\n✗ Validation failed ({check_results['total_errors']} error(s))")
                return 1
            else:
                print("\n✓ All checks passed successfully")
                return 0

        else:
            # Namelist mode — original behavior
            datm_in_path = Path(args.datm_in).resolve()
            if not datm_in_path.exists():
                sys.exit(f"ERROR: datm_in '{datm_in_path}' not found.")

            summary = check_member_namelist(
                datm_in_path, fix_negative=args.fix_negative, verbose=args.verbose,
            )

            print("\n" + "=" * 70)
            print("Summary")
            print("=" * 70)
            if summary["stream_errors"]:
                print(f"Stream file errors: {summary['stream_errors']}")
            if summary["total_errors"]:
                print(f"Forcing file errors: {summary['total_errors']}")
            if not summary["ok"]:
                print("\n✗ Validation failed")
                return 1
            else:
                print("\n✓ All checks passed successfully")
                return 0

    # -----------------------------------------------------------------------
    # Ensemble mode
    # -----------------------------------------------------------------------
    all_results = []

    for i in members:
        if dir_mode:
            base_dir = Path(args.dir).resolve()
            member_dir = base_dir / f"real_{i:05d}"
            print(f"\n{'='*70}")
            print(f"Member: real_{i:05d}")
            print("=" * 70)

            if not member_dir.exists():
                print(f"  ✗ Directory not found: {member_dir}", file=sys.stderr)
                all_results.append({
                    "member_id": f"real_{i:05d}",
                    "ok": False,
                    "total_errors": 1,
                    "stream_errors": 0,
                    "fsds_min": None,
                })
                continue

            check_results = check_member_directory(
                member_dir, variables=variables,
                fix_negative=args.fix_negative, verbose=args.verbose,
            )
            print_forcing_check_results(check_results, verbose=args.verbose)

            fsds_mins = [
                fr["value_checks"]["FSDS"]["min"]
                for fr in check_results["file_results"].values()
                if "FSDS" in fr.get("value_checks", {})
            ]
            all_results.append({
                "member_id": f"real_{i:05d}",
                "ok": check_results["total_errors"] == 0,
                "total_errors": check_results["total_errors"],
                "stream_errors": 0,
                "fsds_min": min(fsds_mins) if fsds_mins else None,
            })

        else:
            # Namelist ensemble mode
            datm_in_path = Path(args.datm_in).resolve()
            member_nml = datm_in_path.parent / f"datm_in_{i:04d}"
            print(f"\n{'='*70}")
            print(f"Member: datm_in_{i:04d}")
            print("=" * 70)

            if not member_nml.exists():
                print(f"  ✗ Not found: {member_nml}", file=sys.stderr)
                all_results.append({
                    "member_id": f"datm_in_{i:04d}",
                    "ok": False,
                    "total_errors": 0,
                    "stream_errors": 1,
                    "fsds_min": None,
                })
                continue

            summary = check_member_namelist(
                member_nml, fix_negative=args.fix_negative, verbose=args.verbose,
            )
            summary["member_id"] = f"datm_in_{i:04d}"
            all_results.append(summary)

    print_ensemble_summary(all_results)

    n_failed = sum(1 for r in all_results if not r["ok"])
    return 0 if n_failed == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
