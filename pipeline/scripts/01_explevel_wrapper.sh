#!/bin/bash
#SBATCH -p short
#SBATCH -t 01:30:00
#SBATCH -c 5
#SBATCH --mem=8G
#SBATCH --array=0-21
#SBATCH -o batch_logs/step01.%j.out 
#SBATCH -e batch_logs/step01.%j.err

# Check if SLURM_ARRAY_TASK_ID is set, if not, define it using the counter variable
if [ -z "$SLURM_ARRAY_TASK_ID" ]; then
    SLURM_ARRAY_TASK_ID=$1
fi 

# CONFIGURE  ENVIRONMENT

source ~/miniconda3/etc/profile.d/conda.sh
conda activate r_env 
source config_pxn.sh # Get global variables from config file

# MAIN 

# Inputs 

metadata=${DSETDIR}/${DSNAME}/metadata_${DSNAME}.csv
tissue_array=($(cut -f1 -d "," ${metadata} | sed 's/\"//g')) # build array of tissues from metadata file
tissue=${tissue_array[${SLURM_ARRAY_TASK_ID}]} # get tissue name form array using job array id
outdir=${OUTDIR}/${GSNAME} # Output directory for results

# Run R script 

echo "---------------------------- Running step 01 - experiment level ----------------------------"
Rscript pxn_01_explevel.R ${DSETDIR} ${GSETDIR} ${outdir} ${DSNAME} ${GSNAME} ${tissue} ${CORES} ${PCOR_OPTION} ${MIN_SAMPLES}
echo "---------------------------- Completed step 01 - experiment level --------------------------"
