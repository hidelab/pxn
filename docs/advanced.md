---
title: Advanced
nav_order: 4
---

# Advanced Uses
{: .no_toc }

## Table of contents
{: .no_toc .text-delta }

1. TOC
{:toc}
---

## Augment mode 

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

4. Run the pipeline as described in the `Getting Started` section, starting from Step 1 - running ./00_prepset-wrapper.sh


## Adding a new background

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

