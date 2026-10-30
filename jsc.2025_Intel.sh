#!/usr/bin/env bash

# Load the JSC module stack when running on JSC. Container images use the
# system GNU/Python/NetCDF tools installed by the Dockerfile instead.
if type module >/dev/null 2>&1; then
    module --force purge
    module use $OTHERSTAGES
    module load Stages/2025
    module load Intel/2024.2.0
    module load ParaStationMPI/5.10.0-1
    module load ESMF/8.7.0
    #
    module load Hypre/2.31.0-cpu
    # Silo not available for the Intel toolchain in Stages/2025 (only gpsmpi/gompi)
    module load Tcl/8.6.14
    #
    module load ecCodes/2.39.0
    module load netCDF/4.9.2
    module load netCDF-Fortran/4.6.1
    module load PnetCDF/1.13.0
    module load cURL/8.7.1
    module load Perl/5.38.2
    module load Python/3.12.3
    module load NCO/5.3.0
    module load CDO/2.4.4
    module load CMake/3.29.3
    module load UnZip/6.0
    #
    module load xarray/2024.11.0
    module load dask/2024.9.1-serial
    module li

    # Export LIB- and INC-NETCDF
    export LIB_NETCDF=${EBROOTNETCDFMINFORTRAN}/lib
    export INC_NETCDF=${EBROOTNETCDFMINFORTRAN}/include
else
    echo "No environment-modules command found; using container toolchain"
    export LIB_NETCDF="${LIB_NETCDF:-/usr/lib}"
    export INC_NETCDF="${INC_NETCDF:-/usr/include}"
fi
