#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Create CORDEX-CMIP6 atmospheric forcing files for eCLM.
"""

import os
import glob
import xarray as xr
import numpy as np
import pandas as pd
import dask
import warnings

# user configuration
root = "/mnt/CORDEX_CMIP6_tmp/sim_data/CORDEX-CMIP6/DD/EUR-12/CLMcom-Hereon/ERA5/evaluation/r1i1p1f1/ICON-CLM-202407-1-1/v1-r1/1hr"
outdir = "/mnt/CORDEX_CMIP6_tmp/user_tmp/spoll/atm-forc/"
variables = ["tas", "rlds", "sfcWind", "pr", "rsds", "psl", "hurs"]
version = "v20240920"
os.makedirs(outdir, exist_ok=True)

rename_map = {
    "tas": "TBOT",
    "rlds": "FLDS",
    "sfcWind": "WIND",
    "rsds": "FSDS",
    "pr": "PRECTmms",
    "psl": "PSRF",
    "hurs": "QBOT"
}

floor_vars = ["rlds", "rsds", "pr"]  # variables to floor instead of round, check https://github.com/WCRP-CORDEX/cordex-cmip6-cmor-tables/blob/main/Tables/CORDEX-CMIP6_1hr.json for time: mean

# dask and warnings
dask.config.set({"array.slicing.split_large_chunks": False})
warnings.filterwarnings("ignore", ".*can't open attribute.*")

###
# helper functions
###
def safe_open_dataset(path):
    try:
        ds = xr.open_dataset(
            path,
            engine="netcdf4",
            decode_cf=False,
            mask_and_scale=False,
            decode_times=False,
            chunks={"time": 720},
        )
        return ds
    except Exception as e:
        print(f"Skipping {path}: {e}")
        return None

def convert_time(ds, varname):
    base_date = np.datetime64("1950-01-01T00:00:00", "ns")
    time_days = ds["time"].values.astype("float64")
    time_days = np.nan_to_num(time_days, nan=0.0)

    if varname in floor_vars:
        hours = np.floor(time_days * 24).astype(int)
    else:
        hours = np.round(time_days * 24).astype(int)

    time_dt = base_date + hours.astype("timedelta64[h]")

    # Drop duplicates
    _, unique_idx = np.unique(time_dt, return_index=True)
    ds = ds.isel(time=unique_idx)
    ds = ds.assign_coords(time=time_dt[unique_idx])
    return ds

def open_and_process(varname):
#    pattern = os.path.join(root, varname, version, f"{varname}_*.nc")
    pattern = os.path.join(root, varname, version, f"{varname}_*196*.nc")
    files = sorted(glob.glob(pattern))
    if not files:
        print(f"No files found for {varname}")
        return None

    # Open with Dask
    datasets = [safe_open_dataset(f) for f in files]
    datasets = [ds for ds in datasets if ds is not None]
    if not datasets:
        return None

    ds_var = xr.combine_by_coords(datasets, combine_attrs="override")
    ds_var = convert_time(ds_var, varname)
    return ds_var

###
# begin of script
###

# load variables
datasets_processed = {}
for var in variables:
    ds_var = open_and_process(var)
    if ds_var is not None:
        datasets_processed[var] = ds_var

if not datasets_processed:
    raise RuntimeError("No valid datasets found.")

# reference dataset
ref_var = "tas" if "tas" in datasets_processed else list(datasets_processed.keys())[0]
ref_ds = datasets_processed[ref_var]

# coordinate variables
coord_vars = {}
for coord_name in ["lat", "lon", "vertices_lon", "vertices_lat", "crs"]:
    if coord_name in ref_ds:
        da = ref_ds[coord_name]
        # If lat/lon carry a time dimension in the source, take first timestep and drop 'time'
        if "time" in da.dims:
            da_no_time = da.isel(time=0, drop=True)
        else:
            da_no_time = da
        coord_vars[coord_name] = da_no_time

###
# monthly loop
###
all_times = pd.to_datetime(ref_ds.time.values)
months = sorted(set(all_times.strftime("%Y%m")))
print("Months found:", months)

for month in months:
    print(f"Processing month {month} ...")
    mask_month = all_times.strftime("%Y%m") == month
    month_start = pd.Timestamp(f"{month[:4]}-{month[4:6]}-01 00:00:00")
    month_vars = {}

    # Process variables individually (otherwise risk of oom) 
    for varname, ds_var in datasets_processed.items():
        ds_month = ds_var.isel(time=np.where(mask_month)[0])

        if ds_month.time.size == 0:
            continue

        # WIND interpolation if needed
        if varname == "sfcWind":
            # Interpolate from 10 m to 2 m, keep roughness height constant (assumption)
            z1, z2, z0 = 10.0, 2.0, 0.01
            ds_month[varname] = ds_month[varname] * (np.log(z2 / z0) / np.log(z1 / z0))

        # ZBOT added later
        month_vars[rename_map[varname]] = ds_month[varname]

    # combine dataset (variables)
    ds_month_all = xr.Dataset(month_vars)

    # add ZBOT
    if "TBOT" in ds_month_all:
        ds_month_all["ZBOT"] = xr.full_like(ds_month_all["TBOT"], 2, dtype=np.float32)

    # add coordinates (needed for interpolation)
    for cname, cdata in coord_vars.items():
        ds_month_all[cname] = cdata
    ds_month_all = ds_month_all.set_coords(["lat", "lon"])

    # set time with metadata (needed for CDO)
    new_time = (pd.to_datetime(ds_month_all.time.values) - month_start) / np.timedelta64(1, "h")
    ds_month_all = ds_month_all.assign_coords(time=("time", new_time))
    ds_month_all["time"].attrs = {
        "standard_name": "time",
        "units": f"hours since {month_start.strftime('%Y-%m-%d %H:%M:%S')}",
        "calendar": "proleptic_gregorian",
        "axis": "T"
    }

    # set global metadata
    ds_month_all.attrs = {
        "author": "Stefan POLL (FZJ)",
        "created_by": "create_atmforcing_CORDEX.py",
        "creation_date": pd.Timestamp.now().strftime("%Y-%m-%d %H:%M:%S"),
    }

    #  write monthly file
    outfile = os.path.join(outdir, f"{month[:4]}-{month[4:6]}.nc")
    print(f"Writing {outfile} ({ds_month_all.time.size} timesteps)...")
    ds_month_all.to_netcdf(
        outfile,
        engine="netcdf4",
        format="NETCDF4",
        unlimited_dims=["time"],
        compute=True,
        encoding={v: {"zlib": True, "complevel": 4} for v in ds_month_all.data_vars}
    )

print("All monthly files written successfully.")
