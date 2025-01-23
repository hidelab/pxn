# PxN System generation pipeline

The system generation consists of three major steps: 

1. Prepating the required input files and objects
2. Calculating pathway correlations at the experiment level
3. Aggregating the pathway estimates across experiments into a combined correlation network

Before jumping into step 1 of the pipeline, it is necessary to generate a standard gene set table, and then modify the config file `config-pdxn.sh` so that all the steps in the pipeline can find the necessary input files and parameters. The following sections describe the process of running the system generation pipeline in detail, starting with modifying the config file.

## Setting up the environment

The PxN pipeline is a collection of R scripts. The required libraries and installation instructions are described in the [Installation Guide](https://github.com/hidelab/PDxN_2.0/blob/main/analysis/pipeline_pdxn_2.0/scripts/example_prep_notebooks/InstallationGuide.ipynb)

## Building the system for single gene set 

### Prepare a standard gene set table 

To run the pipeline with your own gene set you will need to generate a tab-separated file with two columns: `set_name` and `genes`. The `set_name` column indicates the name of the pathway or gene list, and the `genes` column lists the Entrez IDs of the genes in that pathway. If a pathway has 10 genes, this table would have 10 rows for that pathway, one for every gene, where the pathway name is repeated. You can find several examples inside the folder `scripts/example_notebooks` on how to generate this table from unprocessed inputs from public databases such as [MSigSB_v7](https://github.com/hidelab/PDxN_2.0/blob/main/analysis/pipeline_pdxn_2.0/scripts/example_prep_notebooks/prep_geneset-MSigDB_v7.ipynb) or custom gene sets like [Genedex](https://github.com/hidelab/PDxN_2.0/blob/main/analysis/pipeline_pdxn_2.0/scripts/example_prep_notebooks/prep_geneset-Genedex_2024.ipynb). 

Regardless of the starting point, the standard table should look like this:

| set_name          | genes | 
| :---------------- | :------: | 
| Pathway.KEGG_GLYCOLYSIS_GLUCONEOGENESIS        |   55902   | 
| Pathway.KEGG_GLYCOLYSIS_GLUCONEOGENESIS           |   2645   | 
| Pathway.KEGG_GLYCOLYSIS_GLUCONEOGENESIS    |  5232   | 
| Pathway.KEGG_GLYCOLYSIS_GLUCONEOGENESIS |  5230   | 
| Pathway.KEGG_GLYCOLYSIS_GLUCONEOGENESIS        |   5162   | 

### Edit the config file 

The config file contains all the information about where the pipeline inputs are, the resources that will be allocated and other PDxN specific parameters. Every individual script in the pipeline depends on the contents of this config file to run proerly. The first thing you need to do is modify the variables in the config file. DO NOT CHANGE THE VARIABLE NAMES, just modify the file paths or other values (right hand side of the `=` sign). The in-line documentation specifies what to put in each variable. 

Necessary changes:
1. Indicate the path to the input standard gene set table (described above) in the variable `RAW_GSET`.

```{bash}
RAW_GSET='path/to/your/preprocessed/geneset/table/file'
```

2. Give your gene set a short name using `GSNAMEBASE`. If you simply want to use one of our existing gene sets, see the list of opptions available in the [gene set documentation](https://github.com/hidelab/PDxN_2.0/tree/main/analysis/pipeline_pdxn_2.0/input/gene_sets#gene-sets-description).

```{bash}
GSNAMEBASE='MSigDB_v7'
```

3. Define the output folder for your results using `OUTDIR`. We recommend keeping the `output` folder one level up from this `scripts` folder and creating a sub-directory.
```{bash}
OUTDIR='../output/my_run'
```

4. Select the background data set that you want to use with `DSNAME`. To see the list of available data sets refer to the [background data set documentation](https://github.com/hidelab/PDxN_2.0/tree/main/analysis/pipeline_pdxn_2.0/input/gene_expression#background-gene-expression-datasets-description).
```{bash}
DSNAME='gtextoil_gfilter'
```

6. Configure further parameters according to your specific needs. Use the in-line documentation of the config file to learn more about additional settings. 

### Step 1

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
  
### Step 2

#### Server mode

If running on a server, use the script `job-scheduler.sh` to emulate a slurm jobarray to run the first wrapper that processess independently all the tissue groups in the background reference dataset. Simply modify the number of cores directly in the config file (follow its in-code documentation).    

```
nohup ./job-scheduler.sh &
```

*IMPORTANT NOTE ON RESOURCES*: 

The job scheduler will block a given number of cores to launch jobs one after the other. Each of the sequential jobs utilizes multiple cores internally, this parameter is specified in the `config_pdxn.sh` file. You need to take into consideration both numbers when deciding the number of cores you will asign to the job scheduler. For example, if each sequencial job is set to use 8 cores and the scheduler to use 3, then your code will be using 24 cores at any given time. If you set the sequential jobs to use 25 cores each, and the job scheduler to use 4, you would be taking over 120 cores! Be careful and make sure you use only the resources that are needed. A configuration of 8 cores in the sequencial jobs and 3 on the job scheduler strikes a good balance between performance and time. 

*Monitoring jobs*

The job scheduler emulator will create a time-stamped folder inside `scripts/general_logs` that will contain one log file for every 'job', and a general log file. The individual log files show the messages/errors sent to STDOUT by the source script inside the wrapper. The general log file will show a new line every time a new job is launched. You can use this file to track the progress of the code. Once all the jobs have been launched this file sill print a message indicating that it has submitted all jobs. You should also find as many log files as tissues in the selected background dataset. 

#### Cluster mode

If running on a cluster that uses Slurm, modify the job parameters at the top of the script `01_explevel_wrapper.sh` and submit the job by doing:

```
sbatch 01_explevel_wrapper.sh
```

### Step 3

Once the previous step is completed for all tissues/grousp you can run the second step by doing: 

```
./02_combine_wrapper.sh
```

### Exploring the results folder 

The results of a given PXN run will be stored in a subforlder within the directory set in `$OUTDIR` within the config file. If running the pipeline in its single geneset mode, the subfolder will be named according to the following: `[Geneset_name]__[Background_name]`. If running PXN on augment mode, the folder will be named `[Geneset_name]__[Collection_name]__[Network_type]__[Background_name]`. In either case, the folder will contain the following subdirectories:
   - `pathway_activity_tables` - Contains raw pathway activities for all tissues in the background dataset. Each is a `CSV` file with pathways as rows and samples as columns. 
   - `mean_pcor2_barcode_tables` - This folder has the experiment-level pathway correlation estimates resulting from step 1 of the pipeline. There is one file per tissue, containing the pathway pair correlations summarized accross all samples of the same tissue. 
   - `combined_estimates` - This folder contains only three files, which correspond to the aggregation of all tissues in the background dataset. The file `pathway_estimates.tsv` is the most comprehensive one. The other two contain either the P values or the correlations, as indicated in the name of each file.

#### The combined estimates table 

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

## Building the system with a new background gene expression dataset 

Integrating a new background gene expression into the pipeline is as simple as matching the file structure (described below). [This](https://github.com/hidelab/PDxN_2.0/blob/main/analysis/pipeline_pdxn_2.0/scripts/example_prep_notebooks/prep_background-gtex-subset.ipynb) is an example of how to build a custom background dataset by subsetting tissues from the GTex toil dataset. If using a completely different dataset as input, you can use the example to guide your preprocessing choices. 

The important part is that the directory containing your background dataset follows this structure:

- `../input/gene_expression/[DATASETNAME]`   
    - `metadata_[DATASTNAME].csv` CSV file with two columns: tissue name and number of samples
    - `genes_[DATASETNAME].csv` single-column file with Entrez gene ids of all genes in universe(see below)
    - `subgroups`   
        - `tissue_1`   
            - `tissue_1.csv` expression matrix in CSV format (genes x samples) 
        - `tissue_2`    
            - `tissue_2.csv` expression matrix in CSV format (genes x samples) 
         
**Definition of the gene universe**

The file `genes_[DATASETNAME].csv` contains the list of genes in the gene universe of this dataset. This list contains the genes that have at least `min_cts` in at least `min_sam_cts` samples. By default, these parameters are set by default to 3 counts in at least 1 sample.

## Building the system for more than one gene set - augment mode 

This functionality enables the user to compute a network between with two types of nodes (or pathways): core nodes and augmented nodes. The pathways in the core set come from the user's geneset of interest, the one that was used in the previous sections. Augmented nodes are a series of built-in reference pathways that can be incorporated into the analysis of a custom gene set (i.e the core set) to put those pathways of interest in the context of other known sets of genes, such as canonical pathways from MSigDB or drug-response pathways from LINCS (2020 release). This approach can also be further customized with other commonly used gene sets that are not built-in, the process is described [here](https://github.com/hidelab/PDxN_2.0/tree/main/analysis/pipeline_pdxn_2.0/input/augment_sets/README.md). 

To run the pipeline in augment using built-in sets mode follow these steps:

1. Create a collection file (plain text file) inside `input/augment_sets/collections` named [Collection_Name].txt

Browse the list of built-in gene sets inside `input/augment_sets/tables` and write the names of the desired gene sets in the collection file. Only write the portion of the name file preceding the `_pathway_table.csv` suffix. For example, if the gene set file is called `LINCS2020_central_nervous_system_down_pathway_table.csv`, you would write `LINCS2020_central_nervous_system_down` in your collection file. 

2. Enter the name of your collection file without the extension into the config file:

```
COLLECTION='test'
```

3. Specify the type of network that you want to compute in the config file:

```
export NWRK_TYPE='uni_bipartite'
```

There are three potential networks that the system can generate between the core gene set and a collection of reference gene sets:
    a. unipartite: this is a fully-connected network, meaning that all the potential relationships between all nodes are calculated.    
    b. bipartite: this type of network only allows relationships between different types of nodes.
    c. uni-bipartite (recommended): this mode is a custom network where the bipartite network get expanded with a unipartite network of core nodes only.

IMPORTANT: The unipartite mode is the most computationally expensive of all three, since it calculates all pairwaise relationships between nodes. This mode is not recommended for very large gene sets (>15,000 pathways).

4. Run the pipeline as described above, starting from Step 1 - running ./00_prepset-wrapper.sh

## PDxN release notes

| Version          | Notes | 
| :---------------: | :------- | 
|  1.0    |   Original code from [2018 PCxN publication](https://pubmed.ncbi.nlm.nih.gov/29554099/)  | 
|  2.0    |   Re-fctored code to support the use of new background gene expression data set and new gene sets | 
|  2.1   |   Re factored code to optimize runtime by switching to vectorized operations, removed hard-coded parameters, utilize standardized input structures, removed the use of R objects, and condensed code from 4 to 2 main steps | 
| 2.2 |  Incorporated filtering low expressed genes at the top of the pipeline and creating a gene universe to use in downstream steps. Integrated into the gene set pre-processing step a functionality to generate a reference table of pathway pairs along with their overlap coefficients to be used by the entire pipeline. Used the gene universe to filter the genes in the gene set, and added extra heuristics to discard uninformative gene sets and pathway pairs. Because of these changes, this version process a given gene set for each background dataset separately. Which means that the same gene set may have different compositions between two background gene expression datasets. | 
|  2.3      |   Adds the functionality to augment all the pairwise comparison between pathways in a custom gene set with one-sided relationships between pathways in the custom gene set and pathways a set of reference gene sets. The set of reference gene sets includes canonical pathways and tissue-specific drug pathways.  | 





