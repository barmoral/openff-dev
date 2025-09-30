#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=48:00:00
#SBATCH --partition=blanca-shirts
#SBATCH --qos=blanca-shirts
#SBATCH --account=blanca-shirts
#SBATCH --gres=gpu
#SBATCH --mem=4000m
#SBATCH --job-name=FBP_nacl
#SBATCH --output=slurm_codes/FBP_nacl.%j.log

module purge
module avail
ml anaconda
conda activate polymerist-env5

# CASE DETAILS (need to edit each time)
export RESTRAINT=FBP
export SALT="Na Cl"
export CENTERIONS="Na Cl"
export WATER="TIP3P"
export CONCENTRATIONS="0.1 0.5 1.0"
export REPNUM=5
export FORCEK=4148
export FF_WATER='tip3p.offxml'
export FF_IONS='openff-2.1.0.offxml'

#shouldn't need to change
export CASE=${RESTRAINT}_$(echo ${SALT} | tr -d ' ')_${WATER}
export FF_FILES="$FF_WATER $FF_IONS"

python osmotic_sim_dispatch.py -n ${CASE} -s ${SALT} -ci ${CENTERIONS} -w ${WATER} -rn ${REPNUM} -m ${CONCENTRATIONS} -ff ${FF_FILES} -r ${RESTRAINT} -k ${FORCEK} -du nanometer
echo "$CASE $POSTFIX $PDB_INFILE run"

sacct --format=jobid,jobname,cputime,elapsed