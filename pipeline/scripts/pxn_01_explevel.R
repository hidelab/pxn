############################################################
# Script: pdxn_01_explevel.R
#
# Author: Larisa M. Soto
# Adapted from mir_estimates00.R and mir_estimates01.R from
# Pourya Naderi and Sokratis Kariotis
#
# Description:
# This script combines the first two parts of 
# the previous pipeline and performs the following steps for
# a given tissue:
# - Filter gene lists to include only the genes quantified 
#   in the background
# - Calculate gene expression ranks
# - Estimate pathway correlations 
# This script requires input structures generated by the 
# pre-processing script estimates00_2_prep_pdxn_background.ipynb
#
# Arguments:
# - indir: directory containing input gene set objects and background dataset
# - outdir: directory to store outputs
# - dataset_name: name of the background reference dataset
# - geneset_name: name of the gene set 
# - tissue_select: name of the tissue
# - ncores: number of cores 
# - pcor_choice: choice of partial correlation
# - min_samples: min number of samples per tissue
#
# Usage:
# Rscript estimates01.R ../input ../output HGU133plus2_backproc  MSigDB_v7 blood 20 11 0 10 
############################################################
start <- Sys.time() # get start time
# Libraries
suppressPackageStartupMessages({
    library(svd)
    library(corpcor)
    library(future)
    library(furrr)
    library(metap)
    library(data.table)
    library(dplyr) 
})
source('funcs/pdxn_2.0.R')
source('funcs/misc.R')
set.seed(5)

# Input args 

options(stringsAsFactors = F)
cmd_args <- commandArgs(trailingOnly = T)

# User inputs - script args

indir <- cmd_args[1] # Directory with gene expression data
indir_gset <- cmd_args[2] # Directory with gene sets files
output_dir <- cmd_args[3] # Output directory
dataset_name <- cmd_args[4] # Dataset name (subfolder of indir)
geneset_name <- cmd_args[5] # Geneset name (subfolder of indir_gset)
tissue_select <- cmd_args[6] # Name of the tissue to analyze
ncores <- as.numeric(cmd_args[7]) # Number of cores
pcor_choice <- cmd_args[8] # Partial correlation choice
min_samples <- cmd_args[9] # Min number of samples


###############################################
# Main
###############################################

message("\nINPUT ARGS\n")
message("DATASET NAME = ",dataset_name)
message("GENE SET NAME = ",geneset_name)
message("TISSUE = ",tissue_select)
message("NUMBER OF CORES = ",ncores)
message("P. COR CHOICE = ",pcor_choice)
message("MIN SAMPLES = ",min_samples)
message("\n")
# Start parallel multisession

plan(multisession, workers = ncores)

# Prepare dataset structure for PDxN 

## Inputs
input_dir <- file.path(indir,dataset_name)
gset_file <- file.path(indir_gset,geneset_name,"pathway_list.RDS")
pairs_file <- file.path(indir_gset,geneset_name,"pathway_pairs_stats.csv")
cts_file <- file.path(input_dir,"subgroups",tissue_select,paste0(tissue_select,".csv"))
meta_file <- file.path(input_dir,paste0("metadata_",dataset_name,".csv"))

## Outputs 

create_directory(output_dir)
create_directory(file.path(output_dir,"mean_pcor2_barcode_tables"))
create_directory(file.path(output_dir,"pathway_activity_tables"))

outcors_file <- paste0(output_dir,"/mean_pcor2_barcode_tables/",tissue_select,"_cpad_pathcor.tsv")
outsum_file <- paste0(output_dir,"/pathway_activity_tables/",tissue_select,"_path_summaries.csv")
rep.file <- paste(output_dir,"/report.txt",sep = "")
                        
message("Background data set = ", input_dir)
message("Gene set = ", gset_file)
message("Output directory = ",output_dir)

# Load tissue gex data and gene set
#----------------------------------------------------

gs_lst <- readRDS(gset_file)
tissue_exprs <- read.table(cts_file,sep=',',header = TRUE)
tissue_meta <- fread(meta_file,data.table=FALSE,header = TRUE)
pairs <- fread(pairs_file,data.table = FALSE,header = TRUE)
number_of_pathways <- length(gs_lst)

if(ncol(tissue_exprs) < min_samples){
 stop("Not enough samples, cannot process this tissue")
}else{
    message("Tissue = ",tissue_select,"\n",
            "Samples = ",ncol(tissue_exprs),"\n",
            "Pathways = ",number_of_pathways,"\n",
            "Pathway pairs = ",nrow(pairs))
}

# Calculate pathway summaries 
#----------------------------------------------------

## Expression ranks
message('Calculating expression ranks...')
exprs_rnk <- apply(tissue_exprs,2,rank)

## Pathway summaries for disjoint gene sets
message('Estimating pathway summaries...')
summary_dis_list <-list()
for (gs in 1:length(gs_lst)) {
    nm <- names(gs_lst[gs]) # Pathway name
    genes <- as.character(gs_lst[[nm]])
    # Average gene expression ranks for the pathway in all samples in the background
    summary_dis_list[[nm]] <- GetSummary(dat=exprs_rnk,gs=genes,colMeans) 
}
pwy_summs<-do.call(rbind,summary_dis_list) ## Convert pathway summaries to a matrix of size n_pathways x n_samples

# Calculate correlations between pathway pairs
message('Calculating correlations between pairs of pathways...')
res <- pairs %>%
       select(Pathway.A,Pathway.B) %>%
       future_pmap_dfr(~ProcessPair(pathway_a = ..1, pathway_b = ..2)) %>%
       inner_join(pairs,.,by=c('Pathway.A', 'Pathway.B')) 

## Save outputs 
message("Storing outputs...")
write.table(pwy_summs, file = outsum_file, sep = ",", col.names = T, row.names = T, quote = T) # save experiment-level estimates table
write.table(res, file = outcors_file, sep = "\t", col.names = T, row.names = F, quote = F) # save experiment-level estimates table

# Elapsed time 
message('Finished.')
end <- Sys.time() - start 
print(end) 