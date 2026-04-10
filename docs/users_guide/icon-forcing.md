# create_clm-forcing_icon.py

`{note} This page documents the script used to generate **CLM atmospheric forcing files from ICON model output**.`

## Overview

`create_clm-forcing_icon.py` converts atmospheric output from the ICON
model into **forcing data suitable for the Community Land Model (CLM /
eCLM)**.

The script performs the following steps:

1.  Reads ICON NetCDF files
2.  Extracts required meteorological variables
3.  Computes derived quantities required by CLM
4.  Reformats the data structure
5.  Writes CLM-compatible NetCDF forcing files

The resulting dataset contains atmospheric boundary variables required
by CLM land surface simulations.

------------------------------------------------------------------------

## Script Location

mkforcing_icon/create_clm-forcing_icon.py

------------------------------------------------------------------------

# Requirements

## Python Version

Python ≥ 3.8 is recommended.

## Python Dependencies

Main libraries used:

  Package   Purpose
  --------- ------------------------------------
  xarray    NetCDF data processing
  numpy     numerical calculations
  glob      file discovery
  cftime    time handling

------------------------------------------------------------------------

# Input Data

The script processes **ICON atmospheric output files** stored in NetCDF
format.

Input files are collected automatically using a filename pattern.

Example input path:
/path/to/icon/output/ICON_out_EU-R13B5_inst_DOM01_ML_YYYYMM\*.nc

where `MM` represents the month.

## Required Variables

The following variables must be present in the ICON dataset.

  ICON variable   Description
  --------------- ----------------------------
  temp            atmospheric temperature
  qv              specific humidity
  u_10m           zonal wind at 10 m
  v_10m           meridional wind at 10 m
  sob_s           shortwave radiation (down)
  sou_s           shortwave radiation (up)
  thb_s           longwave radiation (down)
  thu_s           longwave radiation (up)
  tot_prec        total precipitation
  pres_msl        mean sea level pressure
  pres_sfc        surface pressure

------------------------------------------------------------------------

# Output

The script produces **monthly CLM forcing NetCDF files**.

Example output:
YYYY-MM.nc

The output dataset contains the following variables.

  CLM Variable   Description
  -------------- --------------------------------------
  FSDS           surface downward shortwave radiation
  FLDS           surface downward longwave radiation
  WIND           near-surface wind speed
  PRECTmms       precipitation rate
  PSRF           surface pressure
  TBOT           near-surface air temperature
  QBOT           near-surface relative humidity
  ZBOT           reference height

Coordinates:

  Coordinate   Description
  ------------ ---------------------
  time         simulation time
  xc           longitude
  yc           latitude
  ni           grid cell index
  nj           singleton dimension

------------------------------------------------------------------------

# Workflow

The script follows several processing steps.

## 1. Collect ICON Files

All monthly ICON files are collected automatically.

filelist = sorted(glob.glob(os.path.join(input_folder, pattern)))

------------------------------------------------------------------------

## 2. Extract Required Variables

Relevant variables are extracted from the dataset.

------------------------------------------------------------------------

## 3. Compute Derived Variables

Several CLM forcing variables must be computed from ICON output.

### Shortwave Radiation

FSDS = sob_s + sou_s

### Longwave Radiation

FLDS = thb_s + thu_s

### Wind Speed

WIND = sqrt(u10² + v10²)

### Precipitation Rate

ICON precipitation is converted from hourly accumulation to a rate:

PRECTmms = tot_prec / 3600

### Relative Humidity

Relative humidity is calculated from temperature, pressure, and specific
humidity.

Example:

RH = 100 \* e / es

where:

-   `e` is vapor pressure
-   `es` is saturation vapor pressure

------------------------------------------------------------------------

## 4. Reshape Data

CLM forcing files expect dimensions:

(time, nj, ni)

ICON data is therefore expanded:

var.expand_dims({"nj":\[1\]}, axis=1)

------------------------------------------------------------------------

## 5. Time Handling

Time values are converted into hours since simulation start.

Example:

hours since YYYY-MM-01 00:00:00

This ensures CF-compliant time encoding.

------------------------------------------------------------------------

## 6. Grid Handling

Longitude and latitude are extracted from ICON grid coordinates.

clon → xc\
clat → yc

Converted from radians to degrees.

------------------------------------------------------------------------

## 7. Domain Truncation

The dataset is truncated to the last `189976` grid cells. Needs to be adopted for others grids.

out = out.isel(ni=slice(-189976, None))

This matches the CLM target grid.

------------------------------------------------------------------------

## 8. Metadata

CF-conventions metadata are added to the output dataset.

Example:

Conventions = CF-1.6

Variable attributes include:

\_FillValue = 1.e20\
missing_value = 1.e20

------------------------------------------------------------------------

## 9. Write Output

Finally, the processed dataset is saved.

Example output filename:
2018-01.nc

------------------------------------------------------------------------

# Example Execution

Run the script with:

python create_clm-forcing_icon.py

Ensure the `input_folder` variable inside the script points to your ICON
output directory.

Example:
/p/scratch/.../icon/

------------------------------------------------------------------------

# Possible Improvements

Future improvements may include:

-   command line argument support
-   configurable grid truncation
-   support for multiple months in one run

