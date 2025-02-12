---
title: Inputs
parent: Getting Started
---

# Input Requirements

To run the pipeline with your own gene set you will need to generate a tab-separated file with two columns: `set_name` and `genes`. The `set_name` column indicates the name of the pathway or gene list, and the `genes` column lists the Entrez IDs of the genes in that pathway. If a pathway has 10 genes, this table would have 10 rows for that pathway, one for every gene, where the pathway name is repeated. You can find several examples inside the folder `scripts/example_notebooks` on how to generate this table from unprocessed inputs from public databases such as [MSigSB_v7](https://github.com/hidelab/PDxN_2.0/blob/main/analysis/pipeline_pdxn_2.0/scripts/example_prep_notebooks/prep_geneset-MSigDB_v7.ipynb) or custom gene sets like [Genedex](https://github.com/hidelab/PDxN_2.0/blob/main/analysis/pipeline_pdxn_2.0/scripts/example_prep_notebooks/prep_geneset-Genedex_2024.ipynb). 

Regardless of the starting point, the standard table should look like this:

| set_name          | genes | 
| :---------------- | :------: | 
| Pathway.KEGG_GLYCOLYSIS_GLUCONEOGENESIS        |   55902   | 
| Pathway.KEGG_GLYCOLYSIS_GLUCONEOGENESIS           |   2645   | 
| Pathway.KEGG_GLYCOLYSIS_GLUCONEOGENESIS    |  5232   | 
| Pathway.KEGG_GLYCOLYSIS_GLUCONEOGENESIS |  5230   | 
| Pathway.KEGG_GLYCOLYSIS_GLUCONEOGENESIS        |   5162   | 
