#!/usr/bin/env bash
set -euo pipefail
if [[ -z "$1" || -z "$2" || -z "$3" || -z "$4" ]]; then
    echo "Usage: $0 MODE YEAR MONTH DOMAINFILE"
    exit 1
fi

MODE=$1
YEAR=$2
MONTH=$3
DOMAINFILE=$4   # "domain.lnd.DE-RuS_DE-RuS.250926.nc"
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"

source "${SCRIPT_DIR}/jsc.2024_Intel.sh"

VENV_DIR="${SCRIPT_DIR}/pyvenv_eclm_atm_forcing_generator"
if [[ ! -d "$VENV_DIR" ]]; then
    python -m venv "$VENV_DIR"
fi
source "${VENV_DIR}/bin/activate"
pip install "${SCRIPT_DIR}"

mkdir -p "${SCRIPT_DIR}/${YEAR}-${MONTH}"
if [[ "$MODE" == "ERA5" ]]; then
    mkdir -p data
    python mkforcing/download_ERA5_input.py \
        --year $YEAR \
        --month $MONTH \
        --dirout data \
        --request "${SCRIPT_DIR}/mkforcing/custom_request_ERA5.py" \
        # --domainfile "${SCRIPT_DIR}/domain.nc"
    mkforcing/prepare_ERA5_input.sh \
        lrenametime=true lmeteo=false \
        lunzip=true wgtcaf=../wgtdis_era5caf_to_domain.nc \
        griddesfile=../domain_griddef.txt iyear=$YEAR \
        imonth=$MONTH \
        pathdata=../data \
        lwgtdis=true lgriddes=true domainfile="${DOMAINFILE}"

else
    python "${SCRIPT_DIR}/mkforcing/download_ERA5_input.py" \
        --year $YEAR --month ${MONTH} \
        --dirout "${SCRIPT_DIR}/cdsapidwn_SEAS5_const" \
        --request "${SCRIPT_DIR}/mkforcing/custom_request_SEAS5_const.py" # \
        # --domainfile $DOMAINFILE
    python "${SCRIPT_DIR}/mkforcing/download_ERA5_input.py" \
        --year ${YEAR} --month ${MONTH} \
        --dirout "${SCRIPT_DIR}/cdsapidwn_SEAS5_24h" \
        --request "${SCRIPT_DIR}/mkforcing/custom_request_SEAS5_24h.py" #  \
        # --domainfile $DOMAINFILE
    python "${SCRIPT_DIR}/mkforcing/download_ERA5_input.py" \
        --year ${YEAR} --month ${MONTH} \
        --dirout "${SCRIPT_DIR}/cdsapidwn_SEAS5_06h" \
        --request "${SCRIPT_DIR}/mkforcing/custom_request_SEAS5_06h.py" #  \
        # --domainfile $DOMAINFILE
    mkdir -p "${SCRIPT_DIR}/cdsapidwn_SEAS5"
    python "${SCRIPT_DIR}/mkforcing/seas5_daily_to_6hourly.py" \
        --const "${SCRIPT_DIR}/cdsapidwn_SEAS5_const/download_era5_${YEAR}_${MONTH}.nc" \
        --daily "${SCRIPT_DIR}/cdsapidwn_SEAS5_24h/download_era5_${YEAR}_${MONTH}.nc" \
        --hourly "${SCRIPT_DIR}/cdsapidwn_SEAS5_06h/download_era5_${YEAR}_${MONTH}.nc" \
        --output "${SCRIPT_DIR}/cdsapidwn_SEAS5/download_era5_${YEAR}_${MONTH}.nc" \
        --frequency 3 --include-hour-zero
    python "${SCRIPT_DIR}/mkforcing/orography_to_elevation.py" \
        "${SCRIPT_DIR}/cdsapidwn_SEAS5/download_era5_${YEAR}_${MONTH}.nc"
    python "${SCRIPT_DIR}/mkforcing/mslp_to_sp.py" \
        "${SCRIPT_DIR}/cdsapidwn_SEAS5/download_era5_${YEAR}_${MONTH}.nc" \
        --elevation-var elevation
    python "${SCRIPT_DIR}/mkforcing/dewpoint_to_specific_humidity.py" \
        "${SCRIPT_DIR}/cdsapidwn_SEAS5/download_era5_${YEAR}_${MONTH}.nc"
    python "${SCRIPT_DIR}/mkforcing/2m_to_10m_conversion.py" \
        "${SCRIPT_DIR}/cdsapidwn_SEAS5/download_era5_${YEAR}_${MONTH}.nc"

    "${SCRIPT_DIR}/mkforcing/prepare_SEAS5_input.sh" \
        lwgtdis=true lgriddes=true \
        domainfile="${DOMAINFILE}" \
        griddesfile="${SCRIPT_DIR}/domain_griddef.txt" \
        wgtcaf="${SCRIPT_DIR}/wgtdis_era5caf_to_domain.nc" \
        iyear=${YEAR} imonth=${MONTH} \
        pathdata="${SCRIPT_DIR}/cdsapidwn_SEAS5"
fi
