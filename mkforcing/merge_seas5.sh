YEAR=$1
NUM_ENS=${2:-51}

mkdir -p "${YEAR}"
for ((i=1; i<=NUM_ENS; i++)); do
  ens=$(printf "%05d" $i)
  files=( ${YEAR}-*/real_${ens}/${YEAR}-*.nc )
  if [ -e "${files[0]}" ]; then
    num_files=${#files[@]}
    echo "Processing ensemble ${ens}: merging ${num_files} files"
    mkdir -p "${YEAR}/real_${ens}"
    cdo -f nc4c mergetime "${files[@]}" "${YEAR}/real_${ens}/${YEAR}-01.nc"
  else
    echo "Ensemble ${ens}: no files found, skipping"
  fi
done
echo "Merging complete for year ${YEAR} with ${NUM_ENS} ensemble members."
