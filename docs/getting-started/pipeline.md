---
title: Pipeline
parent: Getting Started
nav_order: 3
callouts:
  warning:
    title: Warning
    color: red
---

# Running the Pipeline
{: .no_toc }

## Table of contents
{: .no_toc .text-delta }

1. TOC
{:toc}

---

## Overview

The pipeline consists of three major stages:
- **Pre-processing:** Preparing the required input files and directory structure
- **Step 1:** Estimating pathway correlations in each sample
- **Step 2:** Aggregating pathway correlations across 

## Quickstart 

To demonstrate the use of PxN we will simply run the pipeline using the `gtextoil` background with the pre-computed `genedex\_gtextoil` gene set (both described in the [Inputs](https://hidelab.github.io/pxn/docs/getting-started/inputs/) section). This demonstration skips the pre-processing stage, since the gene set is already pre-computed. The first thing you need to do is modify the variables in the config file. Do not change the variable names, just modify the file paths or other values (right hand side of the `=` sign). 

\1. Configure the input and output paramenters in the config file:
   
```
GSNAMEBASE='genedex' # Gene set
DSNAME='gtextoil_gfilter' # Background dataset
OUTDIR='../output/test_run' # Output folder
```

\2. Customize other PxN run parameters according to your system:

```
# Resources
export CORES=10 # number of cores (for part 1)
export CORES_P2=25 # number of cores (for part 2)
```

\3. From the scripts folder, test your configuration by running the first step of the pipeline for one tissue:

```
 ./01_explevel_wrapper.sh 1
```

The standard output should end like this:

```
Storing outputs...
Finished.
Time difference of 18.45294 secs
```

\4. Launch the job scheduler emulator to run the first step for all tissues. Monitor the jobs by looking at the log file inside `general_logs/01_explevel_wrapper_TIMESTAMP/scheduler.log` to see when all jobs are finished (see example below):

```
nohup ./job-sheduler.sh &
```

\5. Once all the jobs are done, run the step 2 of the pipeline. Verify that the code finished successfully by looking at the log file `general_logs/02_combine_wrapper_TIMESTAMP/run.log` (see example below):

```
./02_combine_wrapper.sh
```

## Exploring the results 

The results of a given PXN run will be stored in a subforlder within the directory set in `$OUTDIR` within the config file. If running the pipeline in its single geneset mode, the subfolder will be named according to the following: `[Geneset_name]__[Background_name]`. If running PXN on augment mode, the folder will be named `[Geneset_name]__[Network_type]__[Background_name]`. In either case, the folder will contain the following subdirectories:
   - `pathway_activity_tables` - Contains raw pathway activities for all tissues in the background dataset. Each is a `CSV` file with pathways as rows and samples as columns. 
   - `mean_pcor2_barcode_tables` - This folder has the experiment-level pathway correlation estimates resulting from step 1 of the pipeline. There is one file per tissue, containing the pathway pair correlations summarized accross all samples of the same tissue. 
   - `combined_estimates` - This folder contains only three files, which correspond to the aggregation of all tissues in the background dataset. The file `pathway_estimates.tsv` is the most comprehensive one. The other two contain either the P values or the correlations, as indicated in the name of each file.

**The combined estimates table**

The ultimate output of the PXN pipeline is a table with the pathway correlation estimates that resulted from the aggregation of all tissues in the background dataset. The file looks like this:

| Pathway.A | Pathway.B | p.value | p.adj | PathCor | Overlap.Coefficient | Jacquard.Index | Size.A | Size.B |
| :------------ | :------------- | :-------: | :-------: | :-------: | :-------: | :-------: | :-------: | :-------: |
|AD_Immune_response|Resilience_AD|1|1|-0.078|0|0|42|149|
|AD_Immune_response|AD_CTR|1|1|-0.10|0|0|42|320|
|Resilience_AD|AD_CTR|7.8e-126|4.5e-125|0.26|0.37|0.13|149|320|
|AD_Immune_response|AD_CTR_Astrocytes_down|0.99|1|0.01|0|0|42|306
|Resilience_AD|AD_CTR_Astrocytes_down|2e-138|1.8e-137|0.30|0.03|0.01|149|306|
|AD_CTR|AD_CTR_Astrocytes_down|8.5e-263|1.9e-261|0.41|0.02|0.01|320|306|
|AD_Immune_response|AD_CTR_Astrocytes_up|3.3e-51|8.9e-51|-0.12|0.04|0.05|42|327|
|Resilience_AD|AD_CTR_Astrocytes_up|1.5e-152|1.1e-151|0.31|0.05|0.01|149|327|
|AD_CTR|AD_CTR_Astrocytes_up|1.1e-115|5.9e-115|0.25|0.08|0.04|320|327|


## Running with a custom gene set 

To run PxN with a user-supplied gene set you'll need to:

1. Generate a standard table following the format explained in the Inputs section. Make sure to name your file as `\[GENESETNAME\]_pathway_table.csv`.
2. Place your standard table in the `std_gene_tables` folder.
3. Edit the config file to change the name of the geneset to the name of your geneset (the \[GENESETNAME\] part of the file name you created in step 1).

```
GSNAMEBASE='my-geneset' # Example of GENESETNAME gene set
```

4. Run the input-generation step from the `scripts` folder:

```
./00_prepset_wrapper.sh 
```



This wrapper parses the config file and executes the script `pxn_00_prepset.R` which will process the standard table to meet the pipeline requirements. This script will:

1. Filter based pathways based on their size: Discards pathways that have more than `MAX_GENES` genes or less than `MIN_GENES` genes.
2. Filters the gene set so that it only contains genes that are present in the gene universe of the corresponding backgorund dataset. 
3. Generates the objects required by PDxN inside the `../input` folder:   
    - `pathway_pairs_stats.csv`: a tablle with summary statistics of all pathway pairs passing the Jacquard index threshold (`MAX_JACQ/100`).
    - `pathway_list.RDS`: R object with a list of lists containing the lists of genes in each pathway.
    - `pathway_table.csv`: Human-readable table with the filtered set of gene sets contained in the R object.

You can customize the filtering parameters in the config file. 
  
## Important notes

If running on a server, use the script `job-scheduler.sh` to emulate a slurm jobarray to run the wrapper for step 1 that processess independently all the tissue groups in the background reference dataset. Simply modify the number of cores directly in the config file (follow its in-code documentation).    

**Resource allocation**:

{: .warning }
> The job scheduler will block a given number of cores to launch jobs one after the other. Each of the jobs uses multiple cores internally. You need to take into consideration both numbers when deciding the number of cores you will asign to the job scheduler. For example, if each sequencial job is set to use 8 cores and the scheduler to use 3, then your code will be using 24 cores at any given time. If you set the sequential jobs to use 25 cores each, and the job scheduler to use 4, you would be taking over 120 cores! Be careful and make sure you use only the resources that are needed. A configuration of 8 cores in the sequencial jobs and 3 on the job scheduler strikes a good balance between performance and time. 

**Monitoring jobs**

The job scheduler emulator will create a time-stamped folder inside `scripts/general_logs` that will contain one log file for every 'job', and a general log file. The individual log files show the messages/errors sent to STDOUT by the source script inside the wrapper. The general log file will show a new line every time a new job is launched. You can use this file to track the progress of the code. Once all the jobs have been launched this file sill print a message indicating that it has submitted all jobs. You should also find as many log files as tissues in the selected background dataset. Similarly, the second wrapper creates a log file with the standard output of the R script. 

_Example of general log file_

```
Starting job scheduler...
Max number of cores: 3
Number of jobs to be submitted: 22
[2025-02-16 16:49:34] Launching job 0 of 22 ...
[2025-02-16 16:49:35] Launching job 1 of 22 ...
[2025-02-16 16:49:37] Launching job 2 of 22 ...
[2025-02-16 16:49:48] Launching job 3 of 22 ...
[2025-02-16 16:49:59] Launching job 4 of 22 ...
[2025-02-16 16:50:00] Launching job 5 of 22 ...
[2025-02-16 16:50:11] Launching job 6 of 22 ...
[2025-02-16 16:50:22] Launching job 7 of 22 ...
[2025-02-16 16:50:23] Launching job 8 of 22 ...
[2025-02-16 16:50:34] Launching job 9 of 22 ...
[2025-02-16 16:50:45] Launching job 10 of 22 ...
[2025-02-16 16:50:56] Launching job 11 of 22 ...
[2025-02-16 16:50:57] Launching job 12 of 22 ...
[2025-02-16 16:51:09] Launching job 13 of 22 ...
[2025-02-16 16:51:20] Launching job 14 of 22 ...
[2025-02-16 16:51:21] Launching job 15 of 22 ...
[2025-02-16 16:51:32] Launching job 16 of 22 ...
[2025-02-16 16:51:43] Launching job 17 of 22 ...
[2025-02-16 16:51:44] Launching job 18 of 22 ...
[2025-02-16 16:51:55] Launching job 19 of 22 ...
[2025-02-16 16:52:06] Launching job 20 of 22 ...
[2025-02-16 16:52:07] Launching job 21 of 22 ...
Completed all jobs
```

_Example of step 2 wrapper log file_

```
GENE SET FILE =  ../input/gene_sets/genedex__gtextoil_iBrain/pathway_list.RDS
METADATA FILE =  ../input/gene_expression/gtextoil_iBrain/metadata_gtextoil_iBrain.csv
CORS DIRECTORY =  ../output/test_run/genedex__gtextoil_iBrain/mean_pcor2_barcode_tables
OUTPUT DIRECTORY =  ../output/test_run/genedex__gtextoil_iBrain/combined_estimates
NUMBER OF CORES =  25
Loading experiment-level estimates...
    Number of experiments:  22
Extracting unique set of pathway statistics...
    Number of pathways pairs:  16834
Combining correlation estimates...
Combining P values...
Finished.
Time difference of 3.921468 secs
```
