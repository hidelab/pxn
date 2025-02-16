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

This functionality enables the user to compute a network between with two types of nodes (i.e. gene sets or pathways): _core_ nodes and _augment_ nodes. The nodes in the _core_ set come from the input gene sets file, the one defined at the beginning of the pipeline. The nodes in the _augment_ group come from a built-in reference gene sets file (called _augment_ set)). Classifying the nodes into two categories allows us to derive different kinds of relationships between them. There are three potential networks that the system can generate between _core_ and _augment_ sets:
    a. unipartite: this is a fully-connected network, meaning that all the potential relationships between all nodes are calculated.    
    b. bipartite: this type of network only allows relationships between different types of nodes.
    c. uni-bipartite (recommended): a bipartite network expanded with a unipartite network of core nodes only.

IMPORTANT: The unipartite mode is the most computationally expensive of all three, since it calculates all pairwaise relationships between nodes. This mode is not recommended for very large gene sets (>15,000 pathways).

The pipeline incorporates single or multiple _augment_ set(s) into the standard analysis as a _collection_. Which is simply a text file containing in each line the name of the _augment_ set(s).  

### Inputs 

To showcase this functionality we will use the `test` collection. The `test` collection source file can be found inside the `input` folder at `augment_sets/collections/test.txt`. The file has only 4 lines and looks like this:

```
LINCS2020_autonomic-ganglia_down
LINCS2020_autonomic-ganglia_up
LINCS2020_placenta_down
LINCS2020_placenta_up
```

Each line refers to a different _augment_ set. Their source files are inside `augment_sets/tables`, and each looks like a standard gene set table (described in the [Inputs](https://hidelab.github.io/pxn/docs/getting-started/inputs/#using-custom-gene-sets) section). You do not need to interact with any of these files, they are just shown for clarity. 

### Quickstart

1. Open the config file in a text editor. 

2. Enter the name of your collection file without the extension into the variable `COLLECTION`:

```
COLLECTION='test'
```

3. Specify the type of network type that you want to compute in the variable `NWRK_TYPE`:

```
NWRK_TYPE='uni_bipartite'
```

4. Save your changes to the config file and follow the instructions to run the PxN pipeline with a custom gene set.

### Modifications 

#### Creating a new collection 

Create a collection file (plain text file) inside `pipeline/input/augment_sets/collections` named [CollectionName].txt

Browse the list of built-in gene sets inside `pipeline/input/augment_sets/tables` and write the names of the desired `augment` sets in the collection file (one _augment_ set per line). Only write the portion of the name file preceding the `_pathway_table.csv` suffix. For example, if the gene set file is called `LINCS2020_central_nervous_system_down_pathway_table.csv`, you would write `LINCS2020_central_nervous_system_down` in your collection file. 

#### Adding a new _augment_ set

To add a new _augment_ set to a collection, simply generate a standard table following the same instructions as [using custom gene sets](https://hidelab.github.io/pxn/docs/getting-started/inputs/#using-custom-gene-sets) and place it inside the `tables` folder. Then follow the same procedure described above to add the new _augment_ set to your colection. Make sure to name the file as \[GENESETNAME\]\_pathway\_table.csv, otherwise the pipeline won't recognize it.


## Using a custom background

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

