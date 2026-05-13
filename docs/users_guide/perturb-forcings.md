(perturbforcings)=
# Perturbing Atmospheric Forcings for Ensemble Data Assimilation

```{note}
This step is only needed for eCLM-PDAF data assimilation experiments that
require an ensemble of perturbed atmospheric forcings.
```

The script `mkperturb/perturb_forcings.py` generates an ensemble of perturbed
atmospheric forcing files from a set of monthly NetCDF files (format:
`YYYY-MM.nc`). One perturbed copy per ensemble member is written into
`real_NNNNN/` subdirectories.

## Background

Perturbations follow **Han et al. [2014]**
(https://doi.org/10.1002/2013WR014586), Table 1, and are drawn from a
multivariate normal distribution with cross-correlations between all
four variables applied via Cholesky decomposition:

| Variable            | CLM name   | Noise type                  | Std. dev. |
|---------------------|------------|-----------------------------|-----------|
| Precipitation       | `PRECTmms` | Multiplicative (log-normal) | 0.5       |
| Shortwave radiation | `FSDS`     | Multiplicative (log-normal) | 0.3       |
| Longwave radiation  | `FLDS`     | Additive (normal)           | 20 W/m²   |
| Air temperature     | `TBOT`     | Additive (normal)           | 1 K       |

Cross-correlations between variables (P=Precipitation, SW=Shortwave,
LW=Longwave, AT=Air temperature):

|    |  P   |  SW  |  LW  |  AT  |
|----|------|------|------|------|
| P  |  1.0 | −0.8 |  0.5 |  0.0 |
| SW | −0.8 |  1.0 | −0.5 |  0.4 |
| LW |  0.5 | −0.5 |  1.0 |  0.4 |
| AT |  0.0 |  0.4 |  0.4 |  1.0 |

All other variables in the forcing files (e.g. `PSRF`, `WIND`, `QBOT`,
`ZBOT`) are copied unchanged.

## Prerequisites

The script requires `numpy` and `netCDF4`. The input forcing files must
already be prepared (e.g. by `prepare_ERA5_input.sh`) and named `YYYY-MM.nc`.

## Usage

```bash
python perturb_forcings.py --years <years> --months <months> \
    --num-ensemble <N> --fdir <input_dir> --outdir <output_dir>
```

Years and months accept single values, comma-separated lists, ranges, or
combinations:

```bash
# Single year and month
python perturb_forcings.py --years 2009 --months 1 --num-ensemble 96

# Multiple years and months
python perturb_forcings.py --years 2009,2010,2011 --months 1,2,3

# Range notation
python perturb_forcings.py --years 2009-2011 --months 1-3

# Combined ranges and lists, custom paths
python perturb_forcings.py --years 2009-2011 --months 1-3,7,10-12 \
    --fdir /data/forcings/ --outdir /scratch/ensemble/
```

## Output

For each ensemble member `i` (1-based), a subdirectory `real_NNNNN/` is
created under `--outdir`, containing perturbed copies of all requested
monthly files:

```
<outdir>/
  real_00001/
    2009-01.nc
    2009-02.nc
    ...
  real_00002/
    2009-01.nc
    ...
  ...
```

Each output file carries provenance global attributes (`perturbed_by`,
`perturbed_on_date`) and retains the original attributes under the prefix
`original_attribute_`.

## RNG state and reproducibility

The random number generator state is saved to a JSON file after each run,
allowing the perturbations to be reproduced or continued exactly.

**Default behaviour (`--force-seed`, default):** the RNG is seeded with
`--seed` (default: 42) at the start of each run, overwriting any existing
state file. Runs with the same seed produce identical perturbations.

**Continuing from a previous state (`--no-force-seed`):** the RNG state is
read from the state file saved by a previous run, so the new perturbations
continue the same random sequence rather than restarting it. Useful when
adding more months or ensemble members to an existing set.

```bash
# First run — seeds RNG, saves initial and final state
python perturb_forcings.py --years 2009 --months 1-6 --num-ensemble 96

# Second run — continues from final state of first run
python perturb_forcings.py --years 2009 --months 7-12 --no-force-seed
```

The state file paths default to:
- Initial: `perturb_forcings_rnd_state_initial.json`
- Final: `perturb_forcings_rnd_state_final.json`

A custom path can be set with `--rng-state-file`.

## Perturbing a subset of variables

To perturb only selected variables (e.g. skip longwave radiation):

```bash
python perturb_forcings.py --years 2009-2011 --months 1-12 \
    --variables PRECTmms FSDS TBOT
```

Variables not listed are copied unchanged from the source file.

## All options

| Option                             | Default                                   | Description                                      |
|------------------------------------|-------------------------------------------|--------------------------------------------------|
| `--years`                          | `2009`                                    | Years to process                                 |
| `--months`                         | `1`                                       | Months to process                                |
| `--num-ensemble`                   | `96`                                      | Number of ensemble members                       |
| `--fdir`                           | `./forcings/`                             | Directory containing original forcing files      |
| `--outdir`                         | `./forcings/`                             | Output directory for ensemble subfolders         |
| `--seed`                           | `42`                                      | RNG seed (used when `--force-seed` is active)    |
| `--rng-state-file`                 | `perturb_forcings_rnd_state_initial.json` | Path to RNG state file                           |
| `--variables`                      | all four                                  | Variables to perturb (`PRECTmms FSDS FLDS TBOT`) |
| `--force-seed` / `--no-force-seed` | `--force-seed`                            | Re-seed RNG or continue from saved state         |
