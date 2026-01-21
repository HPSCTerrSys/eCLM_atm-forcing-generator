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

    # Copy netCDF file to a template for remapping (preserves original with all ensemble members)
    cp ${pathdata}/download_era5_${year}_${month}.nc ${tmpdir}/download_era5_${year}_${month}_4rmp.nc

    # Extract the first ensemble member for the template
    # The original file with all 51 ensemble members remains in pathdata for later processing
    ncks --overwrite -d number,0 -O ${tmpdir}/download_era5_${year}_${month}_4rmp.nc ${tmpdir}/download_era5_${year}_${month}_4rmp.nc
    # Remove number and forecast_reference_time dimensions
    ncwa --overwrite -a forecast_reference_time ${tmpdir}/download_era5_${year}_${month}_4rmp.nc ${tmpdir}/download_era5_${year}_${month}_4rmp.nc
    ncwa --overwrite -a number ${tmpdir}/download_era5_${year}_${month}_4rmp.nc ${tmpdir}/download_era5_${year}_${month}_4rmp.nc

    # 1) Renaming variable 'valid_time' to 'time' in $file
    # 2) Renaming dimension 'forecast_period' to 'valid_time' in $file
    # Background: With this naming scheme the CDO remap command below will generate a time variable
    # called "time" in "seconds since 1970-01-01" as for ERA5.
    ncrename -v valid_time,time ${tmpdir}/download_era5_${year}_${month}_4rmp.nc
    ncrename -d forecast_period,valid_time ${tmpdir}/download_era5_${year}_${month}_4rmp.nc

    if $lwgtdis; then
      cdo gendis,${domainfile} ${tmpdir}/download_era5_${year}_${month}_4rmp.nc ${wgtcaf}
    fi

    if $lgriddes; then
      cdo griddes ${domainfile} > ${griddesfile}
    fi

    # TODO: Look into remapping! Now skipped here, but done in the loop below
    # cdo -P ${ompthd} remap,${griddesfile},${wgtcaf} ${tmpdir}/download_era5_${year}_${month}_4rmp.nc ${tmpdir}/rmp_era5_${year}_${month}.nc
  fi

  if $lmerge; then

    # Loop over all 51 ensemble members (indices 0-50)
    for ens in $(seq 0 50); do
      # Format ensemble number as 5-digit with leading zeros (1-based: ens+1)
      ens_num=$(printf "%05d" $((ens + 1)))
      ens_dir=real_${ens_num}
      mkdir -pv ${ens_dir}

      # Copy netCDF file for this ensemble member
      cp ${pathdata}/download_era5_${year}_${month}.nc ${tmpdir}/download_era5_${year}_${month}_ens${ens}.nc

      # Extract the specific ensemble member
      ncks --overwrite -d number,${ens} -O ${tmpdir}/download_era5_${year}_${month}_ens${ens}.nc ${tmpdir}/download_era5_${year}_${month}_ens${ens}.nc
      # Remove number and forecast_reference_time dimensions
      ncwa --overwrite -a forecast_reference_time ${tmpdir}/download_era5_${year}_${month}_ens${ens}.nc ${tmpdir}/download_era5_${year}_${month}_ens${ens}.nc
      ncwa --overwrite -a number ${tmpdir}/download_era5_${year}_${month}_ens${ens}.nc ${tmpdir}/download_era5_${year}_${month}_ens${ens}.nc

      # Rename dimensions/variables for CDO compatibility
      ncrename -v valid_time,time ${tmpdir}/download_era5_${year}_${month}_ens${ens}.nc
      ncrename -d forecast_period,valid_time ${tmpdir}/download_era5_${year}_${month}_ens${ens}.nc

      # Remap this ensemble member
      cdo -P ${ompthd} remap,${griddesfile},${wgtcaf} ${tmpdir}/download_era5_${year}_${month}_ens${ens}.nc ${tmpdir}/rmp_era5_${year}_${month}_ens${ens}.nc

      cdo -P ${ompthd} expr,'WIND=sqrt(u10^2+v10^2)' ${tmpdir}/rmp_era5_${year}_${month}_ens${ens}.nc ${tmpdir}/${year}_${month}_temp.nc # Calculate WIND from u10 and v10
      cdo -f nc4c const,10,${tmpdir}/rmp_era5_${year}_${month}_ens${ens}.nc ${tmpdir}/${year}_${month}_const.nc
      ncpdq -U ${tmpdir}/rmp_era5_${year}_${month}_ens${ens}.nc ${tmpdir}/${year}_${month}_temp2.nc

      cdo merge ${tmpdir}/${year}_${month}_const.nc ${tmpdir}/${year}_${month}_temp2.nc \
          ${tmpdir}/${year}_${month}_temp.nc ${tmpdir}/${year}_${month}_temp4.nc

      ncks -Q -C -x -v hyai,hyam,hybi,hybm ${tmpdir}/${year}_${month}_temp4.nc ${tmpdir}/${year}_${month}_temp5.nc

      # Copy to ensemble-specific directory
      cp ${tmpdir}/${year}_${month}_temp5.nc ${ens_dir}/${year}-${month}.nc

      # Rename variables
      ncrename -v sp,PSRF -v fsds,FSDS -v flds,FLDS -v avg_tprate,PRECTmms -v const,ZBOT -v t10m,TBOT -v q10m,QBOT ${ens_dir}/${year}-${month}.nc

#      ncap2 -O -s 'where(FSDS<0.) FSDS=0' ${ens_dir}/${year}-${month}.nc
      ncatted -O -a units,ZBOT,m,c,"m" ${ens_dir}/${year}-${month}.nc

      ncks -Q -O -h --glb author="${author}" ${ens_dir}/${year}-${month}.nc ${ens_dir}/${year}-${month}.nc
      ncks -Q -O -h --glb contact="${email}" ${ens_dir}/${year}-${month}.nc ${ens_dir}/${year}-${month}.nc

      # Cleanup temporary files for this ensemble member
      rm ${tmpdir}/download_era5_${year}_${month}_ens${ens}.nc
      rm ${tmpdir}/rmp_era5_${year}_${month}_ens${ens}.nc
      rm ${tmpdir}/${year}_${month}_temp*nc ${tmpdir}/${year}_${month}_const.nc
    done
  fi

done
done

