---
title: Pipeline
parent: Getting Started
nav_order: 3
---

# Running the Pipeline
{: .no_toc }

## Table of contents
{: .no_toc .text-delta }

1. TOC
{:toc}

---

## Step 1

Once you've generated this standard table and edited the config file you are ready to run the input-generation step:

```{bash}
./00_prepset_wrapper.sh 
```

This wrapper parses the config file and executes the script `pdxn_00_prepset.R` which will process the standard table to meet the pipeline requirements. This script will:

1. Filter based pathways based on their size: Discards pathways that have more than `MAX_GENES` genes or less than `MIN_GENES` genes.
2. Filters the gene set so that it only contains genes that are present in the gene universe of the corresponding backgorund dataset. 

3. Generates the objects required by PDxN inside the `../input` folder:   
    - `pathway_pairs_stats.csv`: a tablle with summary statistics of all pathway pairs passing the Jacquard index threshold (`MAX_JACQ/100`).
    - `pathway_list.RDS`: R object with a list of lists containing the lists of genes in each pathway.
    - `pathway_table.csv`: Human-readable table with the filtered set of gene sets contained in the R object.
  
## Step 2

### Server mode

If running on a server, use the script `job-scheduler.sh` to emulate a slurm jobarray to run the first wrapper that processess independently all the tissue groups in the background reference dataset. Simply modify the number of cores directly in the config file (follow its in-code documentation).    

```
nohup ./job-scheduler.sh &
```

*IMPORTANT NOTE ON RESOURCES*: 

The job scheduler will block a given number of cores to launch jobs one after the other. Each of the sequential jobs utilizes multiple cores internally, this parameter is specified in the `config_pdxn.sh` file. You need to take into consideration both numbers when deciding the number of cores you will asign to the job scheduler. For example, if each sequencial job is set to use 8 cores and the scheduler to use 3, then your code will be using 24 cores at any given time. If you set the sequential jobs to use 25 cores each, and the job scheduler to use 4, you would be taking over 120 cores! Be careful and make sure you use only the resources that are needed. A configuration of 8 cores in the sequencial jobs and 3 on the job scheduler strikes a good balance between performance and time. 

*Monitoring jobs*

The job scheduler emulator will create a time-stamped folder inside `scripts/general_logs` that will contain one log file for every 'job', and a general log file. The individual log files show the messages/errors sent to STDOUT by the source script inside the wrapper. The general log file will show a new line every time a new job is launched. You can use this file to track the progress of the code. Once all the jobs have been launched this file sill print a message indicating that it has submitted all jobs. You should also find as many log files as tissues in the selected background dataset. 

### Cluster mode

If running on a cluster that uses Slurm, modify the job parameters at the top of the script `01_explevel_wrapper.sh` and submit the job by doing:

```
sbatch 01_explevel_wrapper.sh
```

## Step 3

Once the previous step is completed for all tissues/grousp you can run the second step by doing: 

```
./02_combine_wrapper.sh
```

## Exploring the results 

The results of a given PXN run will be stored in a subforlder within the directory set in `$OUTDIR` within the config file. If running the pipeline in its single geneset mode, the subfolder will be named according to the following: `[Geneset_name]__[Background_name]`. If running PXN on augment mode, the folder will be named `[Geneset_name]__[Collection_name]__[Network_type]__[Background_name]`. In either case, the folder will contain the following subdirectories:
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
