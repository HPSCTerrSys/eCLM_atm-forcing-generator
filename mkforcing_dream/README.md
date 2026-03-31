# Creation of Forcing Data from ICON-DREAM Reanalysis

> **Status:** Experimental – workflow may change in future versions.

## Overview

This project provides a workflow to generate forcing data for **eCLM** based on **ICON-DREAM reanalysis data** from the DWD OpenData platform.

The process includes:
- Downloading ICON-DREAM data
- Remapping the data using DWD ICON tools
- Post-processing adjustments to ensure compatibility with eCLM

## Workflow

1. **Download data**
   - Obtain ICON-DREAM reanalysis data from the DWD OpenData server.

2. **Remapping**
   - Remap the data using DWD ICON tools.
   - Note: This step is planned to be replaced by a CDO-based workflow in future versions.

3. **Post-processing**
   - Adjust variables to match eCLM requirements.
   - Recalculate incoming solar radiation:
     - In ICON-DREAM: accumulated variable
     - In eCLM: instantaneous variable

## Requirements

- DWD ICON tools
- Python environment
- HPC environment (JSC)

## Setup

Before running the script, load the required environment:

```bash
source <environment_file>
source jsc.2024_Intel.sh
```

## Usage

Run the Python script as follows:

```bash
python3 create_forcing.py   --year 2022   --months 6   --base_dir /p/scratch/detectrea2/meurer1   --icon_template create_init_icon_eclm.sh   --tools_workdir /p/scratch/detectrea2/meurer1/dwd_icon_tools   --ingrid /p/scratch/detectrea2/meurer1/dream_grids/ICON-DREAM-EU_grid.nc   --localgrid /p/scratch/detectrea2/meurer1/simexp_DETECT_EUR-3-iic_DWD-ICONglobe_forecast_r1i1p1_FZJ-ICON2024-07-eCLM0-4-0-ParFlow3-14-0_v1/dta/geo/icon/static/EUR-R13B07_2473796_grid_inclbrz_v1.nc   --account detectrea2
```
*(Replace directories and options accordingly.)*

## Credits

- Original DWD ICON tools script: F. Prill  
- Modifications: S. Poll, P. Meurer  
- Python script:
  - Initial version: Pascal Meurer  
  - Further development: ChatGPT v5.2 (LLM-assisted)

- eCLM adaptation code:
  - Based primarily on work by S. Poll  
  - Alternative (slower) version available from P. Meurer  

## Notes

- The workflow is under active development.
- Migration to CDO is planned.
- Further validation of physical consistency (e.g., radiation variables) is ongoing.
