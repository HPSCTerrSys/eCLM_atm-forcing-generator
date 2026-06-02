#!/usr/bin/env python3

import os
import re
import subprocess
import argparse
import shutil

BASE_URL = "https://opendata.dwd.de/climate_environment/REA/ICON-DREAM-EU/hourly"
VARS = ["T", "WS", "QV", "P", "TOT_PREC", "ASWDIFD_S", "ASWDIR_S"]

HEIGHT_VARS = ["T", "WS", "QV", "P"]

# --------------------------------------------------
# Helper
# --------------------------------------------------
def run(cmd): 
    print(">>", " ".join(cmd))
    subprocess.run(cmd, check=True)

def cleanup(paths):
    for p in paths:
        if os.path.isfile(p):
            print(f" deleting file: {p}")
            os.remove(p)
        elif os.path.isdir(p):
            print(f" deleting directory: {p}")
            shutil.rmtree(p, ignore_errors=True)

# --------------------------------------------------
# DOWNLOAD
# --------------------------------------------------
def download_month(year, month, outdir):
    ym = f"{year}{month:02d}"
    os.makedirs(outdir, exist_ok=True)

    files = {}

    for var in VARS:
        fname = f"ICON-DREAM-EU_{ym}_{var}_hourly.grb"
        url = f"{BASE_URL}/{var}/{fname}"
        path = os.path.join(outdir, fname)

        run(["wget", "-c", url, "-O", path])
        files[var] = path

    return files

# --------------------------------------------------
# HEIGHT SELECTION
# --------------------------------------------------
def process_variable(var, infile, outfile, level=74):
    if var in HEIGHT_VARS:
        print(f"{var}: selecting height {level}")
        try:
            run(["cdo", f"sellevel,{level}", infile, outfile])
        except subprocess.CalledProcessError:
            print(f" Level {level} not found → using interpolation")
            run(["cdo", f"intlevel,{level}", infile, outfile])
    else:
        print(f"{var}: no height → copy")
        run(["cdo", "copy", infile, outfile])

# --------------------------------------------------
# MERGE
# --------------------------------------------------
def merge_all(files, outfile):
    print("Merging final file →", outfile)
    run(["cdo", "merge"] + files + [outfile])

# --------------------------------------------------
# PATCH ICON SCRIPT
# --------------------------------------------------
def patch_icon_script(
    template,
    output,
    year,
    month,
    datadir,
    outdir,
    tools_workdir,
    ingrid,
    localgrid,
    account
):
    with open(template) as f:
        content = f.read()

    def replace_var(text, var, value):
        return re.sub(rf"^{var}=.*$", f"{var}={value}", text, flags=re.MULTILINE)

    def replace_slurm(text, key, value):
        return re.sub(rf"^#SBATCH --{key}=.*$", f"#SBATCH --{key}={value}", text, flags=re.MULTILINE)

    content = replace_slurm(content, "account", account)

    content = replace_var(content, "yy", year)
    content = replace_var(content, "mm", f"{month:02d}")
    content = replace_var(content, "tools_WORKDIR", tools_workdir)
    content = replace_var(content, "ICONTOOLS_DIR", "${tools_WORKDIR}/icontools")
    content = replace_var(content, "INGRID", ingrid)
    content = replace_var(content, "LOCALGRID", localgrid)
    content = replace_var(content, "DATADIR", datadir)
    content = replace_var(content, "OUTDIR", outdir)

    datafile = f"${{DATADIR}}/ICON_{year}{month:02d}.grb"
    content = replace_var(content, "DATAFILELIST", datafile)

    if "SLURM_JOB_ID" not in os.environ:
        content = re.sub(r"^srun=.*$", 'srun=""', content, flags=re.MULTILINE)

    # wichtig: sicherstellen dass Verzeichnis existiert
    os.makedirs(os.path.dirname(output), exist_ok=True)

    print("Writing script to:", os.path.abspath(output))

    with open(output, "w") as f:
        f.write(content)

    os.chmod(output, 0o755)

# --------------------------------------------------
# ICON RUN
# --------------------------------------------------
def run_icon(script):
    if "SLURM_JOB_ID" in os.environ:
        run(["sbatch", script])
    else:
        run(["bash", script])

# --------------------------------------------------
# PROCESSING
# --------------------------------------------------
def reconstruct_icon_fsds(da):
    import xarray as xr

    M = da

    H = xr.where(
        (M.time.dt.hour % 3) == 1,
        M,
        xr.where(
            (M.time.dt.hour % 3) == 2,
            2*M - M.shift(time=1),
            xr.where(
                (M.time.dt.hour % 3) == 0,
                3*M - 2*M.shift(time=1),
                M
            )
        )
    )

    return H.clip(min=0).fillna(0)


def process_icon_dream(path_file, y, m, output_path):
    import xarray as xr
    import numpy as np
    import os

    print(f"Processing: {path_file}")

    ds = xr.open_dataset(path_file)

    if "height" in ds.dims:
        ds = ds.sel(height=74.5)

    ncells = ds.sizes["ncells"]
    lon = np.degrees(ds["clon"].values)
    lat = np.degrees(ds["clat"].values)

    TBOT = ds["t"]
    SHUM = ds["q"]
    FSDS_1 = ds["ASWDIFD_S"] + ds["ASWDIR_S"]
    FSDS = reconstruct_icon_fsds(FSDS_1)
    WIND = ds['ws']
    PRECTmms = ds["tp"] / 3600.0
    PSRF = ds['pres']

    ZBOT = xr.DataArray(10*np.ones(ncells), dims=("ncells"))
    def expand_var(var):
        return var.expand_dims({"nj": [1]}, axis=1)

    ds['FSDS'] = expand_var(FSDS)
    ds['WIND'] = expand_var(WIND)
    ds['PRECTmms'] = expand_var(PRECTmms)
    ds['PSRF'] = expand_var(PSRF)
    ds['TBOT'] = expand_var(TBOT)
    ds['SHUM'] = expand_var(SHUM)
    ds['ZBOT'] = expand_var(ZBOT)

    ds['xc'] = expand_var(ds.clon)
    ds['yc'] = expand_var(ds.clat)
    
    ds = ds.set_coords(["xc", "yc"])
    
    keep_vars = ["FSDS", "WIND", "PRECTmms", "PSRF", "TBOT", "SHUM", "ZBOT"]
    ds = ds[keep_vars]

    ds = ds.drop_vars(['nj', 'height', 'clon', 'clat'])

    ds = ds.rename({"ncells": "ni"})

    ds.time.attrs.pop("units", None)

    ds.time.encoding = {
        "units": "seconds since 1970-01-01",
        "calendar": "proleptic_gregorian",
        "dtype": "float64"
    }
    
    target_ni = 189976
    if ds.sizes["ni"] > target_ni:
        ds = ds.isel(ni=slice(-target_ni, None))
        for coord in ["xc", "yc"]:
            ds[coord] = ds[coord].isel(ni=slice(-target_ni, None))

    out_file = f"{y}-{m:02d}.nc"
    ds.to_netcdf(os.path.join(output_path, out_file))

    print(f"Processed monthly file saved to {out_file}")

    print(f"Saved → {out_file}")

# --------------------------------------------------
# MAIN
# --------------------------------------------------
def main():
    parser = argparse.ArgumentParser()

    parser.add_argument("--year", type=int, required=True)
    parser.add_argument("--months", nargs="+", type=int, required=True)

    parser.add_argument("--base_dir", required=True)
    parser.add_argument("--icon_template", required=True)

    parser.add_argument("--tools_workdir", required=True)
    parser.add_argument("--ingrid", required=True)
    parser.add_argument("--localgrid", required=True)
    parser.add_argument("--account", required=True)

    args = parser.parse_args()

    for m in args.months:

        print(f"\n===== {args.year}-{m:02d} =====")

        raw = f"{args.base_dir}/raw/{m:02d}"
        prep = f"{args.base_dir}/prep/{m:02d}"
        merged_dir = f"{args.base_dir}/merged"
        remap_out = f"{args.base_dir}/remap/{m:02d}"

        os.makedirs(prep, exist_ok=True)
        os.makedirs(merged_dir, exist_ok=True)
        os.makedirs(remap_out, exist_ok=True)

        temp_paths = []

        files = download_month(args.year, m, raw)
        temp_paths.append(raw)

        processed_files = []

        for var, path in files.items():
            out = f"{prep}/{var}.grb"
            process_variable(var, path, out, level=74)
            processed_files.append(out)

        temp_paths.append(prep)

        merged = f"{merged_dir}/ICON_{args.year}{m:02d}.grb"
        merge_all(processed_files, merged)

        script = os.path.join(args.base_dir, f"run_icon_{m:02d}.sh")

        patch_icon_script(
            args.icon_template,
            script,
            args.year,
            m,
            merged_dir,
            remap_out,
            args.tools_workdir,
            args.ingrid,
            args.localgrid,
            args.account
        )

        run_icon(script)

        for f in os.listdir(remap_out):
            if f.endswith(".nc"):
                process_icon_dream(
                    os.path.join(remap_out, f),
                    args.year,
                    m,
                    args.base_dir
                )

        cleanup(temp_paths)

    print("\n done")


if __name__ == "__main__":
    main()
