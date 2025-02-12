---
title: Configurations
parent: Getting Started
nav_order: 2
---

# Pipeline Configurations
{: .no_toc }

## Table of contents
{: .no_toc .text-delta }

1. TOC
{:toc}

---

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
