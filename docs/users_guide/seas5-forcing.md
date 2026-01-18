(seas5forcing)=
# eCLM atmospheric forcing based on SEAS5

Start by sourcing the provided environment file

```
source jsc.2024_Intel.sh
```

## Creation of forcing data from SEAS5

Creation of SEAS5 forcing files is adapted from the creation of ERA5
forcing files.

The folder `mkforcing/` contains the scripts that assist the SEAS5
retrieval.

### Download of SEAS5 data

`download_ERA5_input.py` contains a prepared retrieval for the cdsapi python module.

More about cdsapi can be found in [Download of ERA5
data](era5forcing-download).

Usage: Three separate commands have to be executed, one for constant
variables (orography), one for daily variables (e.g. total
precipitation) and one for 6-hourly variables (e.g. temperature).

For each type of output, a dedicated download directory is created.

A main adaption from ERA5 download is the specification of a
customized CDSAPI request using `--request`.

```bash
	python download_ERA5_input.py --year <year> --month <month> --dirout cdsapidwn_SEAS5_const --request ../custom_request_SEAS5_const.py
	python download_ERA5_input.py --year <year> --month <month> --dirout cdsapidwn_SEAS5_24h   --request ../custom_request_SEAS5_24h.py
	python download_ERA5_input.py --year <year> --month <month> --dirout cdsapidwn_SEAS5_06h   --request ../custom_request_SEAS5_06h.py	
```

**Note:** The wrapper script: `./download_ERA5_input_wrapper.sh` is
currently NOT SUPPORTED for SEAS5 download.


### Preparation of SEAS5 data: all variable to 06h

First, we want to have all variables in 6-hourly interval

- from constant: `z` has to be ported (same values as before)
- from daily: `strd`, `ssrd` (thermal, solar, each accumulated) and
  `tp` (total precipitation), these values would have to be
  distributed. For `tp` the value would be divided by four and
  assigned to the four corresponding 6-hour intervals. For the
  radiations, it would be zero in the first/last interval ("night")
  and half the value in the middle two intervals, then converted
  to flux (W/m²) by dividing by the time interval.

The script outputs radiation directly as flux variables (`flds`, `fsds`)
in W/m², so the separate `accumulated_radiation_to_flux.py` step is
not needed for SEAS5 data.

```bash
python seas5_daily_to_6hourly.py --const cdsapidwn_SEAS5_const/download_era5_2026_01.nc --daily cdsapidwn_SEAS5_24h/download_era5_2026_01.nc --hourly cdsapidwn_SEAS5_06h/download_era5_2026_01.nc --output cdsapidwn_SEAS5/download_era5_2026_01.nc
```

### Preparation of SEAS5 data: correct input variables

Steps for preparing SEAS5 data as eCLM input data

1. Orography to elevation (adds a variable `elevation` to the netCDF
   file)

```bash
python orography_to_elevation.py cdsapidwn_SEAS5/download_era5_2026_01.nc
```

2. Mean sea level pressure to surface pressure

```bash
python mslp_to_sp.py cdsapidwn_SEAS5/download_era5_2026_01.nc --elevation-var elevation
```

3. Humidity computed from dewpoint temperature and surface pressure

```
python dewpoint_to_specific_humidity.py cdsapidwn_SEAS5/download_era5_2026_01.nc
```

4. Temperature and Specific Humidity converted from 2m to 10m

```
python 2m_to_10m_conversion.py cdsapidwn_SEAS5/download_era5_2026_01.nc
```

### Preparation of SEAS5 data: Remapping, Data merging, CLM3.5


Check inputs and replace according to your case.

```
sh prepare_SEAS5_input.sh lwgtdis=true lgriddes=true wgtcaf=../wgtdis_era5caf_to_DE-RuS-$(date +%y%m%d).nc griddesfile=../griddes_DE-RuS_$(date +%y%m%d).txt  iyear=2026 imonth=01 author="Johannes KELLER" email="jo.keller@fz-juelich.de" tmpdir=tmpdir pathdata=../cdsapidwn_SEAS5
```
