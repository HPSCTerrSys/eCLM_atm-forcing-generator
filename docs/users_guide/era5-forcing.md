(era5forcing)=
# eCLM atmospheric forcing based on ERA5

Basis: Mainly CDO commands.

By sourcing the provided environment file

```
source jsc.2024_Intel.sh
```

## Creation of forcing data from ERA5

A possible source of atmospheric forcing for CLM (eCLM, CLM5, CLM3.5) is ERA5. It is safer to extract the lowermost level of temperature, humidity and wind of ERA5 instead of taking mixed 2m-values and 10m values. [This internal issue](https://gitlab.jsc.fz-juelich.de/HPSCTerrSys/tsmp-internal-development-tracking/-/issues/36) provides some details. The `download_ERA5_input.py` can be adapted to download another set of quantities.

The folder `mkforcing/` contains three scripts that assist the ERA5 retrieval.

Note: This worfklow is not fully tested.

### Download of ERA5 data

`download_ERA5_input.py` contains a prepared retrieval for the cdsapi python module.
The script requires that cdsapi is installed with a user specific key (API access token).

More information about the installation and access can be found [here](https://cds.climate.copernicus.eu/how-to-api) or alternatively [here](https://github.com/ecmwf/cdsapi?tab=readme-ov-file#install).

Usage:
Either directly:
`python download_ERA5_input.py --year <year> --month <month> --dirout <output_directory>`
Or using the wrapper script:
`./download_ERA5_input_wrapper.sh`
after changing dates and output directory in the `Settings` section inside this wrapper script.

Non-JSC users should adapt the download script to include temperature, specific humidity and horizontal wind speed.

### Preparation of ERA5 data I: Lowermost model level variables (10m altitude)
`extract_ERA5_meteocloud.sh` prepares ERA5 variables form the
lowermost model level (relies on JSC-local files).

Uses level 137, for more information see
https://confluence.ecmwf.int/display/UDOC/L137+model+level+definitions

`extract_ERA5_meteocloud.sh` uses JSC-local grib-input files in
`/p/data1/slmet/met_data/ecmwf/era5/grib/`. Further information:
`/p/data1/slmet/met_data/ecmwf/README.md`.

`extract_ERA5_meteocloud.sh` provides NetCDF files
`meteocloud_YYYY_MM.nc` with lowermost model level atmospheric
variables. The variables from these NetCDF files are used by
`prepare_ERA5_input.sh` in the following ERA5 preparation step.

Usage:
Running the wrapper job
`sbatch extract_ERA5_meteocloud_wrapper.job`
after adapting `year` and `month` loops according to needed dates.

### Preparation of ERA5 data II: Specific humidity computation

For users, who do not have access to the Meteocloud from the previous
section.

For ERA5, specific humidity can be computed from dewpoint temperature
and surface pressure using

```
python dewpoint_to_specific_humidity.py <era5_filename>
```

### Preparation of ERA5 data III: Remapping, Data merging, CLM3.5
`prepare_ERA5_input.sh` prepares ERA5 as an input by remapping the
ERA5 data, changing names and modifying units.

The script is divided into three parts, which could be handled
separately.

1. Remapping
2. Merging the data
3. CLM3.5 specific data preparation.

If remapping is to be used, the remapping weights for the ERA data as
well as the grid definition file of the target domain should be
created beforehand. The following commands can be used to create the
necessary files:

```
cdo gendis,<eclm_domainfile.nc> <era5caf_yyyy_mm.nc> <wgtdis_era5caf_to_domain.nc>
cdo gendis,<eclm_domainfile.nc> <era5meteo_yyyy_mm.nc> <wgtdis_era5meteo_to_domain.nc>
cdo griddes <eclm_domainfile.nc> > <domain_griddef.txt>
```

- `<era5caf_yyyy_mm.nc>`: `caf` stands for "CdsApi Format" and `<era5caf_yyyy_mm.nc>` can be on of the netCDF-files downloaded from the cdsapi.
- `<wgtdis_era5caf_to_domain.nc>` can be chosen, illustrative example:
  `wgtdis_era5caf_to_eur11u-189976.nc`

Usage: `sh prepare_ERA5_input.sh iyear=<year> imonth=<month>
wgtcaf=<wgtcaf> wgtmeteo=<wgtmeteo> griddesfile=<griddesfile>` More
options are available, see script for details.


