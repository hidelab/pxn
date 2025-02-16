#!/bin/bash
#SBATCH -p short
#SBATCH -t 01:30:00
#SBATCH -c 5
#SBATCH --mem=8G
#SBATCH -o batch_logs/step01.%j.out 
#SBATCH -e batch_logs/step01.%j.err

# CONFIGURE  ENVIRONMENT

source ~/miniconda3/etc/profile.d/conda.sh
conda activate r_env 
source config_pxn.sh # Get global variables from config file

# MAIN 

# Create log directory if it doesn't exist
LOG_DIR="general_logs/00_prepset_wrapper_$(date +'%Y-%m-%d_%H-%M-%S')" # Directory to store log files
mkdir -p "${LOG_DIR}"
    
# Inputs 

gsetdir=${GSETDIR}/${GSNAME}
dsetdir=${DSETDIR}/${DSNAME}

# Run R script 

echo "---------------------------- Running step 00 prep set ----------------------------"
Rscript pxn_00_prepset.R ${STD_GSET_TABLE} ${DSNAME} ${GSNAME} ${gsetdir} ${dsetdir} ${GENEUNIV} ${MAX_GENES} ${MIN_GENES} ${MAX_JACQ} ${CORES} ${ASETDIR} ${COLLECTION} ${NWRK_TYPE} >> ${LOG_DIR}/run.log
echo "---------------------------- Completed step 00 prep set --------------------------"
