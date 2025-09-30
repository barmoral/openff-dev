#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=12:00:00
#SBATCH --partition=blanca-shirts
#SBATCH --qos=blanca-shirts
#SBATCH --account=blanca-shirts
#SBATCH --gres=gpu
#SBATCH --mem=6000m
#SBATCH --job-name=HP_mgcl2
#SBATCH --output=slurm_codes/HP_mgcl2.%j.log

module purge
module avail
ml anaconda
conda activate polymerist-env5

# CASE DETAILS (need to edit each time)
export RESTRAINT=HP
export SALT="Mg Cl"
export CENTERIONS="Mg Cl"
export WATER="TIP3P"
export CONCENTRATIONS="0.5"
export REPNUM=6
export FORCEK=0.68
export FF_WATER='tip3p.offxml'
export FF_IONS='openff-2.1.0.offxml'

#shouldn't need to change
export CASE=${RESTRAINT}_$(echo ${SALT} | tr -d ' ')_${WATER}
export FF_FILES="$FF_WATER $FF_IONS"

python osmotic_sim_dispatch.py -n ${CASE} -s ${SALT} -ci ${CENTERIONS} -w ${WATER} -rn ${REPNUM} -m ${CONCENTRATIONS} -ff ${FF_FILES} -r ${RESTRAINT} -k ${FORCEK} -du nanometer
echo "$CASE $POSTFIX $PDB_INFILE run"

sacct --format=jobid,jobname,cputime,elapsed