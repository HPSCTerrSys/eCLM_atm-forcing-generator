#!/usr/bin/env bash
set -euo pipefail
if [[ -z "$1" || -z "$2" || -z "$3" ]]; then
    echo "Usage: $0 MODE YEAR MONTH"
    exit 1
fi

MODE=$1
YEAR=$2
MONTH=$3
DOMAINFILE="domain.lnd.DE-RuS_DE-RuS.250926.nc"

source jsc.2024_Intel.sh 

mkdir -p ${YEAR}-${MONTH}
if [[ "$MODE" == "ERA5" ]]; then
  echo "$MODE not available"
    # mkdir -p data
    # uv run mkforcing/download_ERA5_input.py \
    #     --year $YEAR \
    #     --month $MONTH \
    #     --dirout data \
    #     --request "${HOME}/mkforcing/custom_request_ERA5.py" \
    #     --domainfile "${HOME}/domain.nc"
    # unzip "data/download_era5_${YEAR}_${MONTH}.zip" -d data/
    # uv run mkforcing/dewpoint_to_specific_humidity.py data/data_stream-oper_stepType-instant.nc
    # uv run mkforcing/2m_to_10m_conversion.py data/data_stream-oper_stepType-instant.nc
    # mkforcing/prepare_ERA5_input.sh \
    #     lrenametime=true lmeteo=false \
    #     lunzip=false wgtcaf=../wgtdis_era5caf_to_domain.nc \
    #     griddesfile=../domain_griddef.txt iyear=$YEAR \
    #     imonth=$MONTH \
    #     pathdata=../data \
    #     lwgtdis=true lgriddes=true domainfile="${HOME}/domain.nc"

else
    mkdir cdsapidwn_SEAS5_const
    python mkforcing/download_ERA5_input.py \
        --year $YEAR --month ${MONTH} \
        --dirout cdsapidwn_SEAS5_const \
        --request "mkforcing/custom_request_SEAS5_const.py" \
        --domainfile $DOMAINFILE
    python mkforcing/download_ERA5_input.py \
        --year ${YEAR} --month ${MONTH} \
        --dirout cdsapidwn_SEAS5_24h \
        --request "mkforcing/custom_request_SEAS5_24h.py" \
        --domainfile $DOMAINFILE
    python mkforcing/download_ERA5_input.py \
        --year ${YEAR} --month ${MONTH} \
        --dirout cdsapidwn_SEAS5_06h \
        --request "mkforcing/custom_request_SEAS5_06h.py" \
        --domainfile $DOMAINFILE
    mkdir -p cdsapidwn_SEAS5
    python mkforcing/seas5_daily_to_6hourly.py \
        --const cdsapidwn_SEAS5_const/download_era5_${YEAR}_${MONTH}.nc \
        --daily cdsapidwn_SEAS5_24h/download_era5_${YEAR}_${MONTH}.nc \
        --hourly cdsapidwn_SEAS5_06h/download_era5_${YEAR}_${MONTH}.nc \
        --output cdsapidwn_SEAS5/download_era5_${YEAR}_${MONTH}.nc \
        --frequency 3 --include-hour-zero
    python mkforcing/orography_to_elevation.py \
        cdsapidwn_SEAS5/download_era5_${YEAR}_${MONTH}.nc
    python mkforcing/mslp_to_sp.py \
        cdsapidwn_SEAS5/download_era5_${YEAR}_${MONTH}.nc \
        --elevation-var elevation
    python mkforcing/dewpoint_to_specific_humidity.py \
        cdsapidwn_SEAS5/download_era5_${YEAR}_${MONTH}.nc
    python mkforcing/2m_to_10m_conversion.py \
        cdsapidwn_SEAS5/download_era5_${YEAR}_${MONTH}.nc
    
    mkforcing/prepare_SEAS5_input.sh \
        lwgtdis=true lgriddes=true \
        domainfile=$DOMAINFILE \
        griddesfile=../domain_griddef.txt \
        wgtcaf=../wgtdis_era5caf_to_domain.nc \
        iyear=${YEAR} imonth=${MONTH} \
        pathdata=../cdsapidwn_SEAS5
fi
