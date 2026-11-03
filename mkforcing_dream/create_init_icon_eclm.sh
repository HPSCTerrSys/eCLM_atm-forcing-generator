#!/bin/bash
# ===== SLURM settings ========================================================
#SBATCH --account=detectrea2
#SBATCH --job-name=create_eclm_forcing
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --threads-per-core=2
#SBATCH --output=LOG.create_eclm_forcing.run.%j.o
#SBATCH --error=LOG.create_eclm_forcing.run.%j.o
#SBATCH --time=24:00:00
#=============================================================================

#------------------------------------------yy-----------------------------------
# Run script template for interpolation of inital data onto 
# limited area grid
# 01/2017 : F. Prill, DWD
# 06/2021: SPo modifications for JUWELS
#
# Usage: Submit this PBS runs script with "qsub"
#        Do not forget to fill in the file name settings below!
#
#-----------------------------------------------------------------------------

# ===== srun options ==========================================================
#srun="srun -l --propagate=STACK --cpu_bind=verbose,cores --distribution=block:block"

#test -r /etc/ksh.kshrc && . /etc/ksh.kshrc
set -x

# OpenMP settings
export OMP_SCHEDULE="static"
export OMP_DYNAMIC="false"
export OMP_NUM_THREADS=8

#-----------------------------------------------------------------------------
# MPI variables
# ----------------------------
mpi_root=/p/software/juwels/stages/2024/software/psmpi/5.9.2-1-GCC-12.3.0
no_of_nodes=${SLURM_JOB_NUM_NODES:-1}
mpi_procs_pernode=6
mpi_total_procs=$((no_of_nodes * mpi_procs_pernode))
srun="srun --kill-on-bad-exit=1 --nodes=${SLURM_JOB_NUM_NODES:-1} --ntasks-per-node=${mpi_procs_pernode}"


module purge
module load Stages/2024  GCC/12.3.0  ParaStationMPI/5.9.2-1
module load netCDF-C++4/4.3.1 netCDF-Fortran/4.6.1 HDF5/1.14.2 ecCodes/2.31.0 Szip/.2.1.1 netCDF/4.9.2

module list

# SETTINGS: DIRECTORIES AND INPUT/OUTPUT FILE NAMES --------------------------

# Basic settings:
 yy=2021
 mm=12 #05 #09
 dd=01 #18
 hh=00 

tools_WORKDIR=/p/scratch/detectrea2/meurer1/dwd_icon_tools

ICONTOOLS_DIR=${tools_WORKDIR}/icontools

INGRID=/p/scratch/detectrea2/meurer1/dream_grids/ICON-DREAM-EU_grid.nc

LOCALGRID=/p/scratch/detectrea2/meurer1/simexp_DETECT_EUR-3-iic_DWD-ICONglobe_forecast_r1i1p1_FZJ-ICON2024-07-eCLM0-4-0-ParFlow3-14-0_v1/dta/geo/icon/static/EUR-R13B07_2473796_grid_inclbrz_v1.nc

DATADIR=/p/scratch/detectrea2/meurer1/eCLM_atm-forcing-generator/mkforcing/dream/merged/
DATAFILELIST=$(find ${DATADIR}/*)

OUTDIR=/p/scratch/detectrea2/meurer1/eCLM_atm-forcing-generator/mkforcing/dream/out/

#-----------------------------------------------------------------------------

BINARY_ICONSUB=iconsub
BINARY_REMAP=iconremap
AUXGRID=auxgrid


#-----------------------------------------------------------------------------
# Remap inital data onto local (limited-area) grid
#-----------------------------------------------------------------------------

mkdir -p ${OUTDIR}
cd ${tools_WORKDIR}


set +x

cat >> NAMELIST_ICONREMAP_FIELDS << EOF_2B
!
! U,V - horizontal velocity components
&input_field_nml
 inputname      = "ws"         
 outputname     = "ws"          
 intp_method    = 3
/
! 2m temperature
&input_field_nml
 inputname      = "t"         
 outputname     = "t"          
 intp_method    = 3
/
! pressure
&input_field_nml
 inputname      = "pres"         
 outputname     = "pres"          
 intp_method    = 3
/
! shortwave radiation
&input_field_nml
 inputname      = "ASWDIFD_S"         
 outputname     = "ASWDIFD_S"          
 intp_method    = 3
/
! shortwave radiation
&input_field_nml
 inputname      = "ASWDIR_S"         
 outputname     = "ASWDIR_S"          
 intp_method    = 3
/
! total precipitation
&input_field_nml
 inputname      = "tp"
 outputname     = "tp"
 intp_method    = 3
/
! specific humidity
&input_field_nml
 inputname      = "q"
 outputname     = "q"
 intp_method    = 3
/
EOF_2B

landmask="lsm" ! "FR_LAND"


#-----------------------------------------------------------------------------
# loop over file list:

echo ${DATAFILELIST}
for datafilename in ${DATAFILELIST} ; do

datafile="${datafilename##*/}"  # get filename without path
outdatafile=${datafile%.*}      # get filename without suffix

cat > NAMELIST_ICONREMAP << EOF_2E
&remap_nml
 in_grid_filename  = '${INGRID}'
 in_filename       = '${DATADIR}/${datafile}'
 in_type           = 2
 out_grid_filename = '${LOCALGRID}'
 out_filename      = '${OUTDIR}/${outdatafile}.nc'
 out_type          = 2
 out_filetype      = 4
 l_have3dbuffer    = .false.
 ncstorage_file    = "ncstorage.tmp"
/
EOF_2E


${srun} ${ICONTOOLS_DIR}/${BINARY_REMAP} -vv \
        --remap_nml NAMELIST_ICONREMAP \
        --input_field_nml NAMELIST_ICONREMAP_FIELDS


done

#-----------------------------------------------------------------------------
# clean-up

rm -f ncstorage.tmp*
rm -f nml.log  NAMELIST_SUB NAMELIST_ICONREMAP NAMELIST_ICONREMAP_FIELDS

#-----------------------------------------------------------------------------
exit
#-----------------------------------------------------------------------------

