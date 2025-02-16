---
title: Configurations
parent: Getting Started
nav_order: 2
---

# Pipeline configuration

The pipeline is controlled by the central configuration file `config_pxn.sh`. This config file contains all the information about where the pipeline inputs are, the resources that will be allocated and other PDxN specific parameters. Every individual script in the pipeline depends on the contents of this config file to run proerly. The first thing you need to do is modify the variables in the config file. 

**DO NOT CHANGE THE VARIABLE NAMES**, just modify the file paths or other values (right hand side of the `=` sign).

## Modifications 

The following table outlines the customizable parameters of the config file:

| Variable      | Description | Example | Available Options |
|--------------|------------|---------|-------------------|
| `GSNAMEBASE` | Base name of the gene set | `genedex` | `genedex`, `MSigDBv7`, `MSigDBv6` |
| `DSNAME` | Name of the background reference dataset | `gtextoil_iBrain` | `gtextoil`, `gtextoil_iBrain`, `HGU133plus2` |
| `OUTDIR_NAME` | Name of directory to store outputs | `test_run` | Custom string |
| `COLLECTION` | Optional, name of collection used to augment the custom gene set (augment mode) | `"test"` | Custom string |
| `NWRK_TYPE` | Optional, type of network to calculate (augment) | `"unipartite"` | `unipartite`, `bipartite`, `uni_bipartite` |
| `CORES` | Number of cores (for part 1) | `10` | Integer |
| `CORES_P2` | Number of cores (for part 2) | `25` | Integer |
| `MAX_GENES` | Largest pathway size allowed | `500` | Integer |
| `MIN_GENES` | Smallest pathway size allowed | `20` | Integer |
| `MAX_JACQ` | Max Jacquard index between two pathways | `85` | Integer |
| `PCOR_OPTION` | Choice of partial correlation | `0` | 0 (yes) or 1 (no) |
| `MIN_SAMPLES` | Minimum number of samples per tissue | `10` | Integer |