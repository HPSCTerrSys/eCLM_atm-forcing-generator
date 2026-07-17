#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
perturb_forcings.py — Generate perturbed ensemble atmospheric forcing files.

Reads monthly NetCDF forcing files (format: ``YYYY-MM.nc``) and writes one
perturbed copy per ensemble member into ``real_NNNNN/`` subdirectories.

Perturbations follow Han et al. [2014] (https://doi.org/10.1002/2013WR014586),
Table 1: multiplicative log-normal noise for precipitation and shortwave
radiation, additive normal noise for longwave radiation and air temperature,
with cross-correlations between all four variables.

Usage
-----
  # Single year and month
  perturb_forcings.py --years 2009 --months 1 --num-ensemble 96

  # Multiple years and months (comma-separated list)
  perturb_forcings.py --years 2009,2010,2011 --months 1,2,3

  # Range notation
  perturb_forcings.py --years 2009-2011 --months 1-3

  # Combined ranges and lists, custom paths
  perturb_forcings.py --years 2009-2011 --months 1-3,7,10-12 \\
      --fdir /data/forcings/ --outdir /scratch/ensemble/

  # Continue from existing RNG state instead of reseeding
  perturb_forcings.py --years 2010 --months 1-12 --no-force-seed

  # Perturb only precipitation, shortwave, and temperature (skip longwave)
  perturb_forcings.py --years 2009-2011 --months 1-12 \\
      --variables PRECTmms FSDS TBOT
"""

import numpy as np
import netCDF4 as nc
import os
import json
import datetime
import pathlib
import argparse


def ln_to_n(sd_ln, mean_ln):
    """
    Convert log-normal distribution parameters to the underlying normal parameters.

    Given a log-normal random variable X with mean ``mean_ln`` and standard
    deviation ``sd_ln``, return the mean and standard deviation of the
    underlying normal distribution Y = ln(X).

    Formula from https://en.wikipedia.org/wiki/Log-normal_distribution.

    Parameters
    ----------
    sd_ln : float
        Standard deviation of the log-normal distribution.
    mean_ln : float
        Mean of the log-normal distribution.

    Returns
    -------
    sd_n : float
        Standard deviation of the underlying normal distribution.
    mean_n : float
        Mean of the underlying normal distribution.
    """
    term = 1.0 + sd_ln * sd_ln / mean_ln / mean_ln
    return (np.sqrt(np.log(term)), np.log(mean_ln / np.sqrt(term)))


def n_to_ln(sd_n, mean_n):
    """
    Convert normal distribution parameters to the resulting log-normal parameters.

    Given a normal random variable Y with mean ``mean_n`` and standard
    deviation ``sd_n``, return the mean and standard deviation of X = exp(Y),
    which follows a log-normal distribution.

    Formula from https://en.wikipedia.org/wiki/Log-normal_distribution.

    Parameters
    ----------
    sd_n : float
        Standard deviation of the normal distribution.
    mean_n : float
        Mean of the normal distribution.

    Returns
    -------
    var_ln : float
        Variance of the resulting log-normal distribution.
    mean_ln : float
        Mean of the resulting log-normal distribution.
    """
    return ((np.exp(sd_n * sd_n) - 1.0) * np.exp(2.0 * mean_n + sd_n * sd_n),
            np.exp(mean_n + sd_n * sd_n / 2.0))


def copy_attr_dim(src, dst, usr=None):
    """
    Copy dimensions and global attributes from a source to a destination NetCDF dataset.

    All global attributes from ``src`` are copied to ``dst`` under the prefix
    ``original_attribute_``. Provenance metadata (``perturbed_by``,
    ``perturbed_on_date``) is added to ``dst``.

    Parameters
    ----------
    src : netCDF4.Dataset
        Source dataset to copy dimensions and attributes from.
    dst : netCDF4.Dataset
        Destination dataset to write dimensions and attributes to.
    usr : str, optional
        Username to record in the ``perturbed_by`` attribute. Defaults to the
        ``USER`` environment variable, or ``"unknown"`` if not set.
    """
    # copy attributes
    for name in src.ncattrs():
        dst.setncattr("original_attribute_" + name, src.getncattr(name))
    # copy dimensions
    for name, dimension in src.dimensions.items():
        dst.createDimension(name, len(dimension))
    # Additional attribute
    if usr is None:
        usr = os.environ.get("USER", "unknown")
    dst.setncattr("perturbed_by", usr)
    dst.setncattr("perturbed_on_date",
                  datetime.datetime.today().strftime("%d.%m.%y"))
    # TODO: More attributes, possibly repository-related


ALL_PERTURBED_VARS = ("PRECTmms", "FSDS", "FLDS", "TBOT")


def perturb_nc_file(rng,
                    year=2006,
                    month=1,
                    iensemble=0,
                    fdir="./forcings/",
                    outdir="./forcings/",
                    variables=ALL_PERTURBED_VARS):
    """
    Perturb a single monthly forcing NetCDF file for one ensemble member.

    Reads ``{fdir}/{year}-{month:02d}.nc`` and writes the perturbed output to
    ``{outdir}/real_{iensemble+1:05d}/{year}-{month:02d}.nc``, creating the
    output subdirectory if needed.

    Correlated perturbations (Han et al. [2014], Table 1) are drawn
    from a multivariate normal distribution using Cholesky
    decomposition.

    Any other variables are copied unchanged from the source file.

    Parameters
    ----------
    rng : numpy.random.Generator
        NumPy random number generator (e.g. from ``numpy.random.default_rng``).
    year : int
        Year in the name of the forcing file.
    month : int
        Month in the name of the forcing file.
    iensemble : int
        Zero-based ensemble member index. Output folder uses 1-based numbering.
    fdir : str
        Path to the directory containing the original forcing files.
        Default: ``./forcings/``.
    outdir : str
        Path to the output directory for perturbed ensemble subfolders.
        Default: ``./forcings/``.
    variables : sequence of str, optional
        Forcing variables to perturb. Must be a subset of
        ``("PRECTmms", "FSDS", "FLDS", "TBOT")``. Variables in this set but
        not listed are copied unchanged. Default: all four variables.
    """
    # Perturbation parameters and correlation matrix from:
    # Han et al. [2014] https://doi.org/10.1002/2013WR014586, Table 1
    # Reichle et al. [2007] https://doi.org/10.1029/2006JD008033
    #
    # Variable            | Noise type                  | Std. dev.
    # --------------------|-----------------------------|----------
    # Precipitation       | Multiplicative (log-normal) | 0.5
    # Shortwave radiation | Multiplicative (log-normal) | 0.3
    # Longwave radiation  | Additive (normal)           | 20 W/m²
    # Air temperature     | Additive (normal)           | 1 K
    #
    sd = [ln_to_n(0.5, 1.0)[0], ln_to_n(0.3, 1.0)[0], 20.0, 1.0]
    mean = [ln_to_n(0.5, 1.0)[1], ln_to_n(0.3, 1.0)[1], 0.0, 0.0]

    # Cross-correlation matrix (P=Precip, SW=Shortwave, LW=Longwave, AT=AirTemp):
    #        P     SW    LW    AT
    #   P  [ 1.0, -0.8,  0.5,  0.0]
    #   SW [-0.8,  1.0, -0.5,  0.4]
    #   LW [ 0.5, -0.5,  1.0,  0.4]
    #   AT [ 0.0,  0.4,  0.4,  1.0]
    correl = np.array([[1.0, -0.8, 0.5, 0.0], [-0.8, 1.0, -0.5, 0.4],
                       [0.5, -0.5, 1.0, 0.4], [0.0, 0.4, 0.4, 1.0]])

    fname = pathlib.Path(fdir) / f"{year}-{str(month).zfill(2)}.nc"
    outname = pathlib.Path(outdir) / f"real_{str(iensemble + 1).zfill(5)}" / f"{year}-{str(month).zfill(2)}.nc"

    # Create subdirectories for forcing perturbation
    outname.parent.mkdir(parents=True, exist_ok=True)

    try:
        src = nc.Dataset(fname, "r")
    except FileNotFoundError:
        raise FileNotFoundError(f"Forcing file not found: {fname}") from None

    with src, nc.Dataset(outname, "w") as dst:

        copy_attr_dim(src, dst)

        dim_time = src.dimensions["time"].size
        dim_lat = src.dimensions["lat"].size
        dim_lon = src.dimensions["lon"].size

        # Prepare perturbations:
        # ----------------------

        # Use built-in multivariate, correlated pseudo-number generator
        # directly with Cholesky decomposition
        rnd = rng.multivariate_normal(np.zeros_like(mean), correl,
                                      dim_time * dim_lat * dim_lon,
                                      method="cholesky")

        # Precipitation and ShortWave are log normal distributed perturbations
        # therefore exponential from normal distributed pseudo-random number.

        perturbations = np.zeros_like(rnd)
        perturbations[:, 0] = np.exp(mean[0] + sd[0] * rnd[:, 0])  # Precipitation (multiplicative)
        perturbations[:, 1] = np.exp(mean[1] + sd[1] * rnd[:, 1])  # Shortwave radiation (multiplicative)
        perturbations[:, 2] = mean[2] + rnd[:, 2] * sd[2]          # Longwave radiation (additive)
        perturbations[:, 3] = mean[3] + rnd[:, 3] * sd[3]          # Air temperature (additive)

        # netCDF variables
        # Copy non-perturbed variables:
        # -----------------------------
        for name, var in src.variables.items():
            if name not in variables:
                dst.createVariable(name, var.datatype, var.dimensions)
                dst[name].setncatts(src[name].__dict__)
                dst[name][:] = src[name][:]

        # Add / multiply perturbations
        #-----------------------------

        # Precipitation (multiplicative)
        if "PRECTmms" in variables:
            prectmms = dst.createVariable("PRECTmms",
                                          datatype=np.float64,
                                          dimensions=(
                                              "time",
                                              "lat",
                                              "lon",
                                          ),
                                          fill_value=-9.e+33)
            prectmms.setncatts({"units": u"mm/s", "missing_value": -9.e+33})
            dst.variables["PRECTmms"][:] = (
                src.variables["PRECTmms"][:] * perturbations[:, 0].reshape(
                    src.variables["PRECTmms"][:, :, :].shape))

        # Shortwave radiation (multiplicative)
        if "FSDS" in variables:
            fsds = dst.createVariable("FSDS",
                                      datatype=np.float64,
                                      dimensions=(
                                          "time",
                                          "lat",
                                          "lon",
                                      ),
                                      fill_value=-9.e+33)
            fsds.setncatts({"missing_value": -9.e+33})
            fsds_perturbed = (src.variables["FSDS"][:] *
                              perturbations[:, 1].reshape(src.variables["FSDS"][:, :, :].shape))
            # Clamp small negative values that arise from
            # floating-point noise in the source data.
            #
            # Values more negative than the threshold indicate a
            # genuine problem in the original forcing file and must
            # not be silently discarded.
            _FSDS_NEG_THRESHOLD = -0.01  # W/m²; below this, noise
                                         # clamping is considered
                                         # unsafe
            large_neg = fsds_perturbed[fsds_perturbed < _FSDS_NEG_THRESHOLD]
            if large_neg.size > 0:
                raise ValueError(
                    f"FSDS in {fname} has {large_neg.size} value(s) below "
                    f"{_FSDS_NEG_THRESHOLD} W/m² after perturbation "
                    f"(min={large_neg.min():.4e}). Check the original forcing file."
                )
            dst.variables["FSDS"][:] = np.maximum(0.0, fsds_perturbed)

        # Longwave radiation (additive)
        if "FLDS" in variables:
            flds = dst.createVariable("FLDS",
                                      datatype=np.float64,
                                      dimensions=(
                                          "time",
                                          "lat",
                                          "lon",
                                      ),
                                      fill_value=-9.e+33)
            flds.setncatts({"missing_value": -9.e+33})
            dst.variables["FLDS"][:] = (
                src.variables["FLDS"][:] +
                perturbations[:, 2].reshape(src.variables["FLDS"][:, :, :].shape))

        # Air temperature (additive)
        if "TBOT" in variables:
            tbot = dst.createVariable("TBOT",
                                      datatype=np.float64,
                                      dimensions=(
                                          "time",
                                          "lat",
                                          "lon",
                                      ),
                                      fill_value=9.96921e+36)
            tbot.setncatts({
                "height": u"2",
                "units": u"K",
                "missing_value": -9.e+33
            })
            dst.variables["TBOT"][:, :, :] = (
                src.variables["TBOT"][:, :, :] +
                perturbations[:, 3].reshape(src.variables["TBOT"][:, :, :].shape))


def parse_range(range_str):
    """
    Parse a range string like "1-3,7,10-12" into a list of integers.

    Parameters
    ----------
    range_str : str
        Range specification string.

    Returns
    -------
    list of int
        List of integers.
    """
    result = []
    for part in range_str.split(','):
        if '-' in part:
            start, end = part.split('-')
            result.extend(range(int(start), int(end) + 1))
        else:
            result.append(int(part))
    return sorted(set(result))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description='Perturb atmospheric forcing files.',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Single year and month
  %(prog)s --years 2009 --months 1 --num-ensemble 96

  # Multiple years and months (comma-separated list)
  %(prog)s --years 2009,2010,2011 --months 1,2,3

  # Range notation
  %(prog)s --years 2009-2011 --months 1-3

  # Combined ranges and lists
  %(prog)s --years 2009-2011 --months 1-3,7,10-12
        """)
    parser.add_argument('--years', type=str, default='2009',
                        help='Years to process. Examples: "2009", "2009-2011", "2009,2010,2011", "2009-2011" (default: 2009)')
    parser.add_argument('--months', type=str, default='1',
                        help='Months to process. Examples: "1", "1-3", "1,2,3", "1-3,7,10-12" (default: 1)')
    parser.add_argument('--num-ensemble', type=int, default=96,
                        help='Number of ensemble members (default: 96)')
    parser.add_argument('--fdir', type=str, default='./forcings/',
                        help='Path to directory containing original forcing files (default: ./forcings/)')
    parser.add_argument('--outdir', type=str, default='./forcings/',
                        help='Path to output directory for perturbed ensemble subfolders (default: ./forcings/)')
    parser.add_argument('--seed', type=int, default=42,
                        help='RNG seed used when --force-seed is active (default: 42)')
    parser.add_argument('--rng-state-file', type=str,
                        default='perturb_forcings_rnd_state_initial.json',
                        help='Path to the RNG state file. The final state is saved '
                             'to the same path with "_initial" replaced by "_final", '
                             'or "_final" appended before the extension if "_initial" '
                             'is not in the filename '
                             '(default: perturb_forcings_rnd_state_initial.json)')
    parser.add_argument('--variables', nargs='+',
                        choices=ALL_PERTURBED_VARS,
                        default=list(ALL_PERTURBED_VARS),
                        help='Forcing variables to perturb. Any subset of '
                             'PRECTmms FSDS FLDS TBOT. Variables not listed '
                             'are copied unchanged. '
                             '(default: all four variables)')
    parser.add_argument('--force-seed', action=argparse.BooleanOptionalAction,
                        default=True,
                        help='Force re-seeding the RNG with --seed, overwriting any existing state file (default: true)')
    args = parser.parse_args()

    years = parse_range(args.years)
    months = parse_range(args.months)
    num_ensemble = args.num_ensemble

    # Derive final state file path from initial state file path
    rng_state_file = pathlib.Path(args.rng_state_file)
    _stem = rng_state_file.stem
    _final_stem = _stem.replace("_initial", "_final") if "_initial" in _stem else _stem + "_final"
    rng_state_file_final = rng_state_file.with_name(_final_stem + rng_state_file.suffix)

    # Either seed random number generator or read existing state
    if not os.path.isfile(rng_state_file) or args.force_seed:
        # New Generator instance with seed
        rng = np.random.default_rng(seed=args.seed)

        # Save State of Bit Generator to json-file
        with open(rng_state_file, "w") as f:
            json.dump(rng.bit_generator.state, f)
    else:
        # New Generator instance with state of Bit Generator read from json file
        rng = np.random.default_rng()
        with open(rng_state_file, "r") as f:
            rng.bit_generator.state = json.load(f)

    for ens in range(num_ensemble):
        for y in years:
            for m in months:
                perturb_nc_file(rng, y, m, ens, args.fdir, args.outdir,
                                args.variables)
        print("Done with ensemble " + str(ens + 1))

    # Save final State of Bit Generator to json-file
    with open(rng_state_file_final, "w") as f:
        json.dump(rng.bit_generator.state, f)
