(checkingforcing)=
# Checking Atmospheric Forcing

## Overview

The `check_forcing.py` script validates eCLM atmospheric forcing files
specified in `datm_in` namelist and stream files. It checks for file
existence, required variables as specified in the stream files, and
validates that radiation and precipitation values are non-negative.

## Requirements

```bash
pip install f90nml lxml netCDF4 numpy
```

## Basic Usage

```bash
# Check forcing files (output: errors only)
python checkforcing/check_forcing.py /path/to/datm_in

# Verbose output with detailed information
python checkforcing/check_forcing.py /path/to/datm_in -v

# Automatically fix negative values
python checkforcing/check_forcing.py /path/to/datm_in --fix-negative
```

## What It Checks

1. **Namelist parsing**: Validates `datm_in` and extracts stream file
   references
2. **Stream file parsing**: Reads XML stream files and extracts
   forcing file paths
3. **File existence**: Verifies all referenced forcing files exist
4. **Variable presence**: Checks that NetCDF files contain expected
   variables based on stream mappings
5. **Non-negative validation**: Ensures radiation (`swdn`, `lwdn`) and
   precipitation (`precn`) variables contain no negative values

## Variable Mapping

Stream files define mappings between NetCDF variables (in files) and
eCLM variables (used by model):

```xml
<variableNames>
  FSDS    swdn     <!-- NetCDF variable FSDS maps to eCLM variable swdn -->
  FLDS    lwdn
  PRECTmms  precn
</variableNames>
```

## Output Modes

### Quiet Mode (default)
Shows only errors and final status:
```
✗ 2022-10.nc:
    FSDS (eCLM: swdn): 12/8760 values are negative (min=-0.0034)

✗ Validation failed
```

### Verbose Mode (`-v`)
Shows complete information including:
- Namelist details
- Stream file contents
- All forcing files checked
- Variable mappings
- Detailed statistics

### Fix Mode (`--fix-negative`)
Automatically corrects negative values:
- Sets negative values to zero
- Updates file history attribute
- Re-validates after fixing
```
✓ Fixed 2022-10.nc: FSDS(12)

✓ All checks passed successfully
```

## Exit Codes

- `0`: All checks passed
- `1`: Validation errors found
