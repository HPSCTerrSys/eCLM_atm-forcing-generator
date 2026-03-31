# Creation of Forcing Data from ICON-DREAM Reanalysis

> ⚠️ **Status:** Experimental – workflow may change in future versions.

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
   - ⚠️ Note: This step is planned to be replaced by a CDO-based workflow in future versions.

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
python <script_name>.py [options]
```

*(Replace `<script_name>` and options accordingly.)*

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
