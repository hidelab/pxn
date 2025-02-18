--- 
title: Inputs
parent: Getting Started
nav_order: 1
---

# Pipeline Inputs
{: .no_toc }

## Table of contents
{: .no_toc .text-delta }

1. TOC
{:toc}

---

The PxN pipeline needs two main components: (1) a gene expression background, and (2) a gene set table. PxN comes with a series of gene expression datasets and gene set files ready to use. It is also possible for the user to provide a custom gene set and/or a custom background dataset. Changing the background dataset requires a several pre-processing steps, and its expalined in the `Advanced` section. The sections below describe the built-in PxN datasets and outline the steps needed to incorporate a custom gene set into the pipeline. 

## Overview

PxN comes with two main background datasets and several gene sets. The pre-processed gene sets come from three major sources:

- MSigDB version 7 C2 collection (internally tagged as `MSigDBv7`)
- MSigDB version 6 C2 collection (internally tagged as `MSigDBv6`)
- Genedex: a manually curated database of Alzheimer\'s Disease facets (internaly tagged as `genedex`)

The background datasets are:

_GTex Toil_

This dataset is part of the [UCSC Toil RNAseq Recompute Compendium](https://xenabrowser.net/datapages/?host=https%3A%2F%2Ftoil.xenahubs.net). Raw files from the `TOIL_RSEM_norm_count` were downloaded from the [Xena browser](https://xenabrowser.net/datapages/?cohort=GTEX&removeHub=https%3A%2F%2Fxena.treehouse.gi.ucsc.edu%3A443). The internal tag for this dataset is `gtextoil`, it includes 19,561 genes across 7,847 samples from 54 different tissues. A subset of this dataset containing only brain and immune tissues is available with the tag `gtextoil_iBrain`.

_Microarray Barcode_

This dataset was obtained from the set of legacy files of the orginal PCxN version. It is provided for reference and comparison with the previous version of the pipeline. The internal tag for this dataset is `HGU133plus2`, it contains 20,590 genes across 3,207 samples from 72 tissues.  This dataset was obtained from the input files of the orginal [PCxN](https://pubmed.ncbi.nlm.nih.gov/29554099/). 

## Data Download 

PxN built-in data folder be downloaded as a `.tgz` file from Zenodo. To make sure that it integrates into the pipeline ecosystem, the file needs to be expanded inside the `input` folder under `pipeline`.   

For example, if you cloned the repo into `~/pxn` your directory strcucture should look like this:

    - ~/pxn/pipeline
    	- input/
    	- scripts/  
    	- output/

1. Download the data from Zenodo [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.14879142.svg)](https://doi.org/10.5281/zenodo.14879142)
2. Move the dowloaded file into the `-/pxn/pipeline/input` folder
3. While inside the `-/pxn/pipeline/input` folder, expand the file:

```
cd ~/pxn/pipeline/input # Enter the input directory
tar -xvzf input.tgz # Expand the tar file
rm input.tgz # Delete the tar file after expanding
```

## Contents

Each subfolder comes with a `README.md` file providing extensive documentation about its contents. In summary, the `input` folder contains the following subdirectories:

| Directory   | Contents  |
|------------|----------|
| augment_sets | Described in the `Advanced` section |
| gene_expression | Processed gene expression data |
| gene_sets | Processed gene set objects |
| std_gene_tables | Reference gene set tables |
    
### The `gene_expression` folder

This folder includes post-processed gene expression datasets structured in the format required for PxN. The processing steps undertaken to generate these datasets are:
- Keeping only genes expressed with at least 3 counts in at least 1 sample (they constitute the _gene universe_)
- Discarding tissues with less than 10 samples.

List of datasets:

- `gtextoil`- This dataset uses as base the GTex toil dataset.
    - Gene universe size: 19561
    - Number of samples: 7847
    - Number of tissues: 54

- `gtextoil_iBrain` - This is a subset of the `gtextoil` dataset that includes only brain and immune tissues. 
    - Gene universe size: 18994
    - Number of samples: 2746
    - Number of tissues: 22

- `HGU133plus`

This dataset uses as base the HGU133plus microarray barcode dataset.
    - Gene universe size: 20590
    - Number of samples: 3207
    - Number of tissues: 72

### The `gene_sets` folder

Each reference gene set gets processed for a particular gene expression background to ensure that only genes expressed in the background dataset are included in the gene set used for analysis. This results in unique background-gene set combinations that are labelled using the gene set and background dataset internal tags. 

List of gene sets:

- `MSigDBv7\__gtextoil`- This geneset contains the pathways from MSigDB version 7 filtered for the gene universe of the GTex toil (gfilter) dataset.
  - Pathways: 1186
  - Pathway pairs: 702689

- `MSigDBv6\__gtextoil`- This geneset contains the pathways from MSigDB version 6 filtered for the gene universe of the GTex toil (gfilter) dataset.
  - Pathways: 851
  - Pathway pairs: 361628

- `genedex\__gtextoil` - This gene set contains the gene lists of AD facets available in Genedex (accessed Nov 2024), processed for the GTex toil (gfilter) dataset.
  - Pathways: 242
  - Pathway pairs: 16834

- `MSigDBv7\__HGU133plus2` - This geneset contains the pathways from MSigDB version 7 filtered for the gene universe of the HGU133plus2 barcode (gfilter) dataset.
  - Pathways: 1199
  - Pathway pairs: 718189

- `MSigDBv6\__gtextoil` - This geneset contains the pathways from MSigDB version 6 filtered for the gene universe of the GTex toil (gfilter) dataset.
  - Pathways: 851
  - Pathway pairs: 361628

### The `std_gene_tables` folder

This directory contains a series of standard tables prepared for different reference gene sets. These files are used as the source to prepare the gene set files in the `gene_sets` folder. There are three tables included by default in the pipeline, one for each source database. The naming convention is `\[GENESETNAME\]\_pathway\_table.csv` This is where custom user-specific standard tables (see below) should be place to ensure consistency.

## Creating custom gene sets

To run the pipeline with your own gene set you will need to generate a standard table, which is basically a tab-separated file with two columns: `set_name` and `genes`. This file should be place inside the `std_gene_tables` folder following the file naming convention of standard tables (described above). The `set_name` column indicates the name of the pathway or gene list, and the `genes` column lists the Entrez IDs of the genes in that pathway. If a pathway has 10 genes, this table would have 10 rows for that pathway, one for every gene, where the pathway name is repeated. You can find several examples inside the folder `scripts/example_notebooks` on how to generate this table from unprocessed inputs from public databases such as [MSigSB](https://github.com/hidelab/pxn/pipeline/scripts/example_prep_notebooks/prep_geneset-MSigDB_v7.ipynb) or custom gene sets like [Genedex](https://github.com/hidelab/pxn/pipeline/scripts/example_prep_notebooks/prep_geneset-Genedex_2024.ipynb). 

Regardless of the starting point, the standard table should look like this:

| set_name          | genes | 
| :---------------- | :------: | 
| Pathway.KEGG_GLYCOLYSIS_GLUCONEOGENESIS        |   55902   | 
| Pathway.KEGG_GLYCOLYSIS_GLUCONEOGENESIS           |   2645   | 
| Pathway.KEGG_GLYCOLYSIS_GLUCONEOGENESIS    |  5232   | 
| Pathway.KEGG_GLYCOLYSIS_GLUCONEOGENESIS |  5230   | 
| Pathway.KEGG_GLYCOLYSIS_GLUCONEOGENESIS        |   5162   | 
