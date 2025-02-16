#!/bin/bash

# =============================================================================
# Script: job-scheduler.sh
# Author: Larisa M. Soto (lmorale3@bidmc.harvard.edu)
# Description:
# This script acts as a job scheduler to run a specified script that 
# that takes as input a sequential run number ID (similar to a jobarray in slurm
# on multiple input files concurrently, utilizing a maximum number of available
# CPU cores up to a user defined limit. The script ensures that at no point 
# exceed the specified number of parallel jobs. Each job's output is logged 
# with timestamps. The script dynamically determines the number of available CPU 
# cores and uses up to the max cores if available. 
#
# Usage:
# 1. Edit the Parameter section to specify:
#    the script to be executed, the log directory, 
# 2. Run the script:
#    nohup ./job-scheduler.sh &
#
# Parameters:
# - SCRIPT: The script to run with each input file.
# - LOG_DIR: The directory to store log files.
# - NUM_JOBS: Number of jobs to run (upper limit of 1-based counter)
#
# Functions:
# - get_running_jobs: Counts the number of currently running jobs.
# - run_job: Runs a job for a given input file and logs the output with timestamps.
#
# =============================================================================

source config_pxn.sh # Get input variables for the run 

########################################
# Parameters
########################################

SCRIPT="01_explevel_wrapper.sh"    # The actual script you want to run
LOG_PARENT_DIR="general_logs"      # Parent directory to store log folders
MAX_CORES=3                        # Max number of cores used by the job scheduler

# To specify the number of cores manually 

#NUM_JOBS=25 

# To build the job scheduler based on the number of input files in a given directory

INPUT_DIR=${DSETDIR}/${DSNAME}/subgroups  # Directory containing the input files for each script run
input_files=("$INPUT_DIR"/*)
NUM_JOBS=${#input_files[@]} # This is the upper limit of the counter variable

########################################
# Setup
########################################

# Get the number of available cores, default to 10 if more cores are available
AVAILABLE_CORES=$(nproc)
MAX_CORES=$(( AVAILABLE_CORES < MAX_CORES ? AVAILABLE_CORES : MAX_CORES ))

# Create log directory if it doesn't exist
LOG_DIR="$LOG_PARENT_DIR/$(basename "$SCRIPT" .sh)_$(date +'%Y-%m-%d_%H-%M-%S')" # Directory to store log files
mkdir -p "${LOG_PARENT_DIR}"
mkdir -p "$LOG_DIR"

# Print startup message
echo "Starting job scheduler..."
echo "Log file directory: $LOG_DIR"
echo "Scheduler log file: ${LOG_DIR}/scheduler.log"

echo "Starting job scheduler..." >> ${LOG_DIR}/scheduler.log
echo "Max number of cores: $MAX_CORES" >> ${LOG_DIR}/scheduler.log
echo "Number of jobs to be submitted: $NUM_JOBS" >> ${LOG_DIR}/scheduler.log

########################################
# Functions
########################################

get_num_log_files() {
  local num_logs=$(find "$LOG_DIR" -maxdepth 1 -type f | wc -l)
  echo "$num_logs"
}

get_running_jobs() {
  # Uses ps and grep to count the number of running instances of the script
  echo $(ps aux | grep "$SCRIPT" | grep -v grep | wc -l)
}

run_job() {
  local i=$1               # The counter variable (index of the input file)
  local log_file=$2        # The log file name
  # Run the script with the input file and redirect output to the log file with timestamp
  {
    echo "[$(date +'%Y-%m-%d %H:%M:%S')] Starting job for sample number $i"
    bash "$SCRIPT" "$i"
    echo "[$(date +'%Y-%m-%d %H:%M:%S')] Completed job for sample number $i"
  } &> "$log_file" &
}

########################################
# Main 
########################################

i=0
while [ $(( $(get_num_log_files) - 1 )) -lt "$NUM_JOBS" ]; do

  # Ensure the number of running jobs does not exceed the maximum allowed cores
  while [ $(get_running_jobs) -ge $MAX_CORES ]; do
    #echo "Waiting to submit next job..." >> ${LOG_DIR}/scheduler.log
    sleep 10 # Wait for 10 seconds before checking again
  done

  # Construct the log file name
  log_file="$LOG_DIR/run_$i.log"
  
  # Run the job with the current index and log file name
  echo "[$(date +'%Y-%m-%d %H:%M:%S')] Launching job $i of $NUM_JOBS ..." >> ${LOG_DIR}/scheduler.log
  run_job "$i" "$log_file"
  
  # Increment the index for the next input file
  ((i++))
  
  # Wait for a second before starting the next job
  sleep 1
done

echo "Completed all jobs"
echo "Completed all jobs" >> ${LOG_DIR}/scheduler.log
