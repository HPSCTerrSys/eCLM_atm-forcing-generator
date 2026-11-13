#!/usr/bin/env bash

# Load the JSC module stack when running on JSC. Container images use the
# system GNU/Python/NetCDF tools installed by the Dockerfile instead.
#
# NOTE: Module version strings below are taken from Stages/2025 as a starting
# point.  Verify each one with `module avail <name>` under Stages/2026 and
# update accordingly before using this file.
if type module >/dev/null 2>&1; then
    module --force purge
    module use $OTHERSTAGES
    module load Stages/2026
    module load Intel
    module load ParaStationMPI
    module load ESMF
    #
    module load Hypre
    module load Tcl
    #
    module load ecCodes
    module load netCDF
    module load netCDF-Fortran
    module load PnetCDF
    module load cURL
    module load Perl
    module load Python
    module load NCO
    module load CDO
    module load CMake
    module load UnZip
    #
    module load xarray
    module load dask
    module li

    # Export LIB- and INC-NETCDF
    export LIB_NETCDF=${EBROOTNETCDFMINFORTRAN}/lib
    export INC_NETCDF=${EBROOTNETCDFMINFORTRAN}/include
else
    echo "No environment-modules command found; using container toolchain"
    export LIB_NETCDF="${LIB_NETCDF:-/usr/lib}"
    export INC_NETCDF="${INC_NETCDF:-/usr/include}"
fi
