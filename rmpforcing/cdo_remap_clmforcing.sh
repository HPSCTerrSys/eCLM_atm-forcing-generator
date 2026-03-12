# Shell script to remap clm5 forcing
# 

# module load
ml Stages/2024  Intel/2023.2.1  ParaStationMPI/5.9.2-1 CDO/2.3.0 NCO

# forcing
clm_forcdir=/p/data1/jibg31/DETECT_newgrid/ERA5_reana/dataOut/3hr
clm_domainfile=/p/scratch/cslts/poll1/sim/euro-cordex/tsmp2_workflow-engine/dta/geo/eclm/static/domain.lnd.ICON-11_maskedhalo.230302.nc
outdir=/p/scratch/cslts/poll1/data/eclm_forc/eur11u/

#
cdo griddes $clm_domainfile > ${outdir}/clmgrid.txt
#clm_forcnam=${clm_forcdir}/2022-07.nc

# 
#cdo -P 2 remapdis,${outdir}/clmgrid.txt ${clm_forcnam} ${outdir}/2022-07.nc
# cdo remap,clmgrid.txt,${outdir}/rmp_weights.nc ${clm_forcnam} ${outdir}/2022-07_ss.nc

# loop over month
date_sta='2007-01-01'
date_end='2017-12-01'
endt=$(date '+%s' -d "$date_end")
idate="$date_sta"

while [[ $(date +%s -d $idate) -le $endt ]]; do
#   ${outdir}/${idate%-*}.nc
   forcname=${idate%-*}.nc
   # remap file
   cdo -P 2 remapdis,${outdir}/clmgrid.txt ${clm_forcdir}/$forcname ${outdir}/$forcname
   idate=$(date '+%Y-%m-%d' -d "$idate +1 month")
done

