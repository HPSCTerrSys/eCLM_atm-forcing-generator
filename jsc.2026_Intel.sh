#!/usr/bin/env bash

# Load the JSC module stack when running on JSC. Container images use the
# system GNU/Python/NetCDF tools installed by the Dockerfile instead.
#
if type module >/dev/null 2>&1; then
    module --force purge
    module use $OTHERSTAGES
    module load Stages/2026
    module load Intel/2025.2.0
    module load ParaStationMPI/5.13.0-1
    module load ESMF/8.9.0
    #
    module load Hypre/3.0.0-cpu
    module load Tcl/9.0.1
    #
    module load ecCodes/2.44.0
    module load netCDF/4.9.3
    module load netCDF-Fortran/4.6.2
    module load PnetCDF/1.14.1
    module load cURL/8.14.1
    module load Perl/5.40.2
    module load Python/3.13.5
    module load NCO/5.3.6
    module load CDO/2.5.4
    module load CMake/4.0.3
    module load UnZip/6.0
    #
    module load xarray/2025.10.0
    module load dask/2025.9.1
    module li

    # Export LIB- and INC-NETCDF
    export LIB_NETCDF=${EBROOTNETCDFMINFORTRAN}/lib
    export INC_NETCDF=${EBROOTNETCDFMINFORTRAN}/include
else
    echo "No environment-modules command found; using container toolchain"
    export LIB_NETCDF="${LIB_NETCDF:-/usr/lib}"
    export INC_NETCDF="${INC_NETCDF:-/usr/include}"
fi
