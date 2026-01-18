#!/usr/bin/env bash
set -eo pipefail

# default values of parameters
lrmp=true
lmerge=true
lwgtdis=false			# Switch for creating wgtdis file in script
lgriddes=false			# Switch for creating griddes file in script
ompthd=1

# TSMP2/eclm
pathdata=./
domainfile=../domain.lnd.DE-RuS_DE-RuS.250926.nc
wgtcaf=/p/scratch/cslts/poll1/sim/euro-cordex/tsmp2_wfe_eur-11u/dta/rmp_gridwgts/wgtdis_era5caf_to_eur11u-189976.nc
# wgtmeteo=/p/scratch/cslts/poll1/sim/euro-cordex/tsmp2_wfe_eur-11u/dta/rmp_gridwgts/wgtdis_era5meteo_to_eur11u-189976.nc
griddesfile=/p/scratch/cslts/poll1/sim/euro-cordex/tsmp2_wfe_eur-11u/dta/rmp_gridwgts/griddes_eur-11u_189976.txt

iyear=2017
imonth=07
tmpdir=tmpdir
wrkdir=""
author="Default AUTHOR"
email="d.fault@fz-juelich.de"

# Function to parse input
parse_arguments() {
    for arg in "$@"; do
        key="${arg%%=*}"
        value="${arg#*=}"

        case "$key" in
            lrmp) lrmp="$value" ;;
            lmerge) lmerge="$value" ;;
            lwgtdis) lwgtdis="$value" ;;
            lgriddes) lgriddes="$value" ;;
            ompthd) ompthd="$value" ;;
            pathdata) pathdata="$value" ;;
            wgtcaf) wgtcaf="$value" ;;
            # wgtmeteo) wgtmeteo="$value" ;;
            griddesfile) griddesfile="$value" ;;
            tmpdir) tmpdir="$value" ;;
            wrkdir) wrkdir="$value" ;;
            imonth) imonth="$value" ;;
            iyear) iyear="$value" ;;
            author) author="$value" ;;
            email) email="$value" ;;
            *) echo "Warning: Unknown parameter: $key" ;;
        esac
    done
}

# Call the function to parse the input arguments
# Users needs to make sure for consistent input
parse_arguments "$@"

#
#cd $wrkdir
#mkdir -pv $tmpdir

#
for year in ${iyear}
do
for month in ${imonth}
do

  # Go into working directory and create temporary directory
  if [ -z ${wrkdir} ];then
    wrkdir=${iyear}-${imonth}
  fi
  cd $wrkdir
  mkdir -pv $tmpdir

  if $lrmp; then

    # Copy netCDF file
    cp ${pathdata}/download_era5_${year}_${month}.nc ${tmpdir}

    # Extract the first ensemble member
    ncks --overwrite -d number,0 -O ${tmpdir}/download_era5_${year}_${month}.nc ${tmpdir}/download_era5_${year}_${month}.nc
    # Remove number and forecast_reference_time dimensions
    ncwa --overwrite -a forecast_reference_time ${tmpdir}/download_era5_${year}_${month}.nc ${tmpdir}/download_era5_${year}_${month}.nc

    # Renaming variable 'valid_time' to 'time' in $file
    ncrename -v valid_time,time ${tmpdir}/download_era5_${year}_${month}.nc

    if $lwgtdis; then
      cdo gendis,${domainfile} ${tmpdir}/download_era5_${year}_${month}.nc ${wgtcaf}
    fi

    if $lgriddes; then
      cdo griddes ${domainfile} > ${griddesfile}
    fi

    cdo -P ${ompthd} remap,${griddesfile},${wgtcaf} ${tmpdir}/download_era5_${year}_${month}.nc ${tmpdir}/rmp_era5_${year}_${month}.nc
  fi

  if $lmerge; then

    cdo -P ${ompthd} expr,'WIND=sqrt(u10^2+v10^2)' ${tmpdir}/rmp_era5_${year}_${month}.nc ${tmpdir}/${year}_${month}_temp.nc # Calculate WIND from u10 and v10
    cdo -f nc4c const,10,${tmpdir}/rmp_era5_${year}_${month}.nc ${tmpdir}/${year}_${month}_const.nc
    ncpdq -U ${tmpdir}/rmp_era5_${year}_${month}.nc ${tmpdir}/${year}_${month}_temp2.nc

    cdo merge ${tmpdir}/${year}_${month}_const.nc ${tmpdir}/${year}_${month}_temp2.nc \
        ${tmpdir}/${year}_${month}_temp.nc ${tmpdir}/${year}_${month}_temp4.nc

    ncks -C -x -v hyai,hyam,hybi,hybm ${tmpdir}/${year}_${month}_temp4.nc ${tmpdir}/${year}_${month}_temp5.nc

    # Simply copy the file
    cp ${tmpdir}/${year}_${month}_temp5.nc  ${year}-${month}.nc

    # Rename variables
    ncrename -v sp,PSRF -v fsds,FSDS -v flds,FLDS -v tp,PRECTmms -v const,ZBOT -v t10m,TBOT -v q10m,QBOT ${year}-${month}.nc

#    ncap2 -O -s 'where(FSDS<0.) FSDS=0' ${year}_${month}.nc
    ncatted -O -a units,ZBOT,m,c,"m" ${year}-${month}.nc

    ncks -O -h --glb author="${author}" ${year}-${month}.nc
    ncks -O -h --glb contact="${email}" ${year}-${month}.nc

    rm ${tmpdir}/${year}_${month}_temp*nc ${tmpdir}/${year}_${month}_const.nc
  fi

done
done

