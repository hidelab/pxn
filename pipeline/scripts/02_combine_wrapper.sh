#!/bin/bash
#SBATCH -p short
#SBATCH -t 03:00:00
#SBATCH --mem=24G
#SBATCH --cpus-per-task=25
#SBATCH -o batch_logs/step02.%j.out 
#SBATCH -e batch_logs/step02.%j.err

# CONFIGURE  ENVIRONMENT

source ~/miniconda3/etc/profile.d/conda.sh
conda activate r_env 
source config_pxn.sh # Get global variables from config file

# MAIN 

# Inputs

geneset=${GSETDIR}/${GSNAME}/pathway_list.RDS # Geneset file
metadata=${DSETDIR}/${DSNAME}/metadata_${DSNAME}.csv # Dataset metadata file
pcordir=${OUTDIR}/${GSNAME}/mean_pcor2_barcode_tables # Directory with partial correlations
outdir=${OUTDIR}/${GSNAME}/combined_estimates # Output directory for results
LOG_PARENT_DIR="general_logs"      # Parent directory to store log folders

# Create log directory if it doesn't exist
LOG_DIR="$LOG_PARENT_DIR/02_combine_wrapper_$(date +'%Y-%m-%d_%H-%M-%S')" # Directory to store log files
mkdir -p "${LOG_PARENT_DIR}"
mkdir -p "${LOG_DIR}"

# Run R script 

echo "---------------------------- Running step 02 - combine ----------------------------"
Rscript pxn_02_combine.R ${geneset} ${metadata} ${pcordir} ${outdir} ${CORES_P2} >> ${LOG_DIR}/run.log
echo "---------------------------- Completed step 02 - combine --------------------------"
