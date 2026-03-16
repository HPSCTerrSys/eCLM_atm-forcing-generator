import xarray as xr
import numpy as np
import glob
import os
import cftime
from datetime import timedelta 

def calc_relative_humidity(temp, pres, qv):
    """
    Calculate relative humidity in % from temperature (K), pressure (Pa), and specific humidity (kg/kg)
    Using the standard formula: RH = 100 * e / es
    """
    # Saturation vapor pressure over water (Pa)
    T_C = temp - 273.15
    es = 610.94 * np.exp((17.625 * T_C) / (T_C + 243.04))
    # Actual vapor pressure
    e = qv * pres / (0.622 + 0.378 * qv)
    RH = 100.0 * e / es
    return RH

# --- Step 1: collect all ICON files for the month ---
month = 1  # 
mm_str = f"{month:02d}"
input_folder = f"/p/scratch/cslts/poll1/sim/paper/wfe_eur-11_revsetup_icon-eclm-pfl/dta/simres/iconeclmparflow_2018{mm_str}01/out/icon/"  # adjust to your path
filelist = sorted(glob.glob(os.path.join(input_folder, f"ICON_out_EU-R13B5_inst_DOM01_ML_2018{mm_str}*.nc")))

if not filelist:
    raise FileNotFoundError("No ICON files found. Check your path and filename pattern.")

# --- Step 2: open all files as one dataset ---
vars_needed = [
    "temp", "qv", "u_10m", "v_10m",
    "sob_s", "sou_s", "thb_s", "thu_s",
    "tot_prec", "pres_msl"
]

ds = xr.open_mfdataset(
    filelist,
    combine="nested",
    concat_dim="time",
    data_vars=vars_needed,
    coords="minimal",
    compat="override"
)

temp = ds["temp"].isel(height=59)
qv = ds["qv"].isel(height=59)

u10 = ds["u_10m"].isel(height_3=0)
v10 = ds["v_10m"].isel(height_3=0)

# --- Step 3: extract lon/lat from ICON ---
lon = np.degrees(ds["clon"].values)
lat = np.degrees(ds["clat"].values)
ncells = ds.dims["ncells"]

# --- Step 4: calculate derived variables ---

# Shortwave downward surface flux
FSDS = ds["sob_s"] + ds["sou_s"]

# Longwave downward surface flux
FLDS = ds["thb_s"] + ds["thu_s"]

# Near-surface wind magnitude
#WIND = np.sqrt(ds["u_10m"]**2 + ds["v_10m"]**2)
WIND = np.sqrt(u10**2 + v10**2)

# Precipitation rate [kg/m²/s] (convert from hourly kg/m²)
PRECTmms = ds["tot_prec"] / 3600.0

# Sea level pressure
PSRF = ds["pres_msl"]

# Temperature
TBOT = temp

# Relative humidity
#QBOT = calc_relative_humidity(ds["temp"], ds["pres_sfc"], ds["qv"])
QBOT = calc_relative_humidity(temp, ds["pres_sfc"], qv)

# ZBOT constant 10
ZBOT = xr.DataArray(10*np.ones((ds.dims["time"], ncells)), dims=("time","ncells"))

# --- Step 5: expand dimensions to match (time, nj=1, ni=ncells) ---
def expand_var(var):
    return var.expand_dims({"nj": [1]}, axis=1)

FSDS = expand_var(FSDS)
FLDS = expand_var(FLDS)
WIND = expand_var(WIND)
PRECTmms = expand_var(PRECTmms)
PSRF = expand_var(PSRF)
TBOT = expand_var(TBOT)
QBOT = expand_var(QBOT)
ZBOT = ZBOT.expand_dims({"nj": [1]}, axis=1)

icon_time0 = ds["time"].values[0]
time_hours = (ds["time"].values - icon_time0) * 24.0
time_hours = np.round(time_hours, 6).astype("float64")  # fractional hours

y = int(str(int(icon_time0))[:4])
m = int(str(int(icon_time0))[4:6])
d = int(str(int(icon_time0))[6:8])
start_time_str = f"hours since {y}-{m:02d}-{d:02d} 00:00:00"

# --- Step 6: create dataset ---
out = xr.Dataset(
    {
        "FSDS": (("time", "nj", "ni"), FSDS.values),
        "FLDS": (("time", "nj", "ni"), FLDS.values),
        "WIND": (("time", "nj", "ni"), WIND.values),
        "PRECTmms": (("time", "nj", "ni"), PRECTmms.values),
        "PSRF": (("time", "nj", "ni"), PSRF.values),
        "TBOT": (("time", "nj", "ni"), TBOT.values),
        "QBOT": (("time", "nj", "ni"), QBOT.values),
        "ZBOT": (("time", "nj", "ni"), ZBOT.values),
    },
    coords={
        "time": ("time", time_hours),
        "xc": (("nj", "ni"), lon[None, :]),
        "yc": (("nj", "ni"), lat[None, :]),
        "nj": [1],
        "ni": np.arange(ncells),
    }
)

# important step
out["time"].attrs.update({
    "standard_name": "time",
    "axis": "T",
    "units": start_time_str,          # copy into attrs for NetCDF
    "calendar": "proleptic_gregorian"
})

out["time"].encoding.clear()
out["time"].encoding.update({
    "units": f"hours since {start_time_str}",
    "calendar": "proleptic_gregorian",
    "dtype": "float64"
})

# --- Step 7: truncate to last 189,976 points for all variables and coordinates ---
target_ni = 189976
if out.dims["ni"] > target_ni:
    out = out.isel(ni=slice(-target_ni, None))
    for coord in ["xc", "yc"]:
        out[coord] = out[coord].isel(ni=slice(-target_ni, None))

# --- Step 8: add metadata (example from processed file) ---
for var in ["FSDS","FLDS","WIND","PRECTmms","PSRF","QBOT","ZBOT"]:
    out[var].attrs["_FillValue"] = 1.e20
    out[var].attrs["missing_value"] = 1.e20

out.attrs["Conventions"] = "CF-1.6"
out.attrs["author"] = "Stefan POLL (FZJ)"
out.attrs["created_by"] = "ICON -> CLM forcing conversion script"

# --- Step 9: save output ---
out_file = f"2018-{mm_str}.nc"
out.to_netcdf(out_file)
print(f"Processed monthly file saved to {out_file}")
