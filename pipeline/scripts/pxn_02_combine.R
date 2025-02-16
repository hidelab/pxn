##########################################################################################################
# Script: pdxn_02_combine.R
#
# Author: Larisa M. Soto
# Adapted from mir_estimates02.R and mir_estimates03.R from
# Pourya Naderi and Sokratis Kariotis
#
# Description:
#
# This script aggregates the experiment-level estimates into aggregated estimates for every pathway pair 
# It weights the correlations and the pvalues by the number of samples in each experiment to calculate
# the pathway-level estimates. A final FDR adjustment is applied to the combined p values. The resulting 
# files consitute the final PCxN network of patwhay-pathway correlations in the backrground dataset. 

# Arguments:
#
# 1. geneset_file: Geneset file
# 2. gsovlp_file: Geneset overlaps file
# 3. meta_file: Dataset metadata file
# 4. cors_dir: Directory with partial correlations
# 5. res_folder: Output directory 
# 6. ncores: Number of cores 
#
# Usage:
#
# Rscript estimates002.R ../input/gene_sets/MSigDB_v7/pathway_list.RDS \
#                        ../input/gene_sets/MSigDB_v7/ pathway_pairs_overlaps.RDS \
#                        ../input/gene_expression/HGU133plus2_backproc/metadata_HGU133plus2_backproc.csv \
#                        ../output/HGU133plus2_backproc_MSigDB_v7/mean_pcor2_barcode_tables \
#                        ../output/HGU133plus2_backproc_MSigDB_v7/combined_estimates
#
##########################################################################################################

globalstart <- Sys.time()
suppressPackageStartupMessages({
    library(parallel)
    library(data.table)
    library(metap)
    library(dplyr) 
    library(tidyr) 
})
source('funcs/pdxn_2.0.R')
source('funcs/misc.R')

# Input settings 

set.seed(5)
options(stringsAsFactors = F)
options(future.globals.maxSize = 4000 * 1024^2)
cmd_args <- commandArgs(trailingOnly = T)

# User inputs - script args

geneset_file <- cmd_args[1]  # Geneset file
meta_file <- cmd_args[2] # Dataset metadata file
cors_dir <- cmd_args[3] # Directory with partial correlations
res_folder <- cmd_args[4] # Output directory for results
ncores <- cmd_args[5] # Number of cores 

cat("GENE SET FILE = ",geneset_file)
cat("\nMETADATA FILE = ",meta_file)
cat("\nCORS DIRECTORY = ",cors_dir)
cat("\nOUTPUT DIRECTORY = ",res_folder)
cat("\nNUMBER OF CORES = ",ncores)

create_directory(res_folder)

########################################################################
# Main
########################################################################

# Load data
gs_lst <- readRDS(geneset_file)
meta <- read.table(meta_file,sep=",",header = T)
cors_files <- list.files(cors_dir,full.names = TRUE)

# Load experiment level estimates
cat("\nLoading experiment-level estimates...")
cat("\n    Number of experiments: ",length(cors_files))

exp_cors <- cors_files %>%
            mclapply(.,function(path){
                df <- data.table::fread(path,header=TRUE) %>%
                      mutate(experiment=sub("_cpad_pathcor.tsv","",basename(path)))
                return(df)                
            },mc.cores=ncores) 
cors_df <- do.call(rbind,exp_cors)

# Calculate weights
meta <- meta %>%
        mutate(weight = nsamples/sum(nsamples))
exp_weights <- meta$weight
names(exp_weights)<-meta$subgroup

# Extract pathway statistics
cat("\nExtracting unique set of pathway statistics...")
res <- cors_df %>%
        select(-estimate,-n,-statistic,-p.value,-experiment) %>%
        distinct()
cat("\n    Number of pathways pairs: ", nrow(res))

# Combine correlation estimates
cat("\nCombining correlation estimates...")
cors_wide <- experiments_to_wide(cors_df,values_col = "estimate") # Re-format data 

if(length(names(exp_weights))!=length(colnames(cors_wide))){
    stop('Missing tissue outputs. Verify that the previous step was completed successfully for all tissues.')
}

path_cor <- as.matrix(cors_wide[,names(exp_weights)]) %*% exp_weights # Calculate weighted average

# Re-format data into a data frame with the two pathways in each pair split into two columns
pcor_df <- data.frame(Pathway.A=gsub(";.*","",rownames(path_cor)), 
                      Pathway.B=gsub(".*;","",rownames(path_cor)),
                      PathCor=path_cor)
rownames(pcor_df)<-NULL
        
# Store partial table
write.table(pcor_df,
            file = file.path(res_folder,"pathway_correlations.tsv"),
            sep="\t",col.names=TRUE,row.names=FALSE,quote=FALSE)
        
rm(exp_cors,path_cor,cors_wide) # Cleanup environment

# Combine P values
cat("\nCombining P values...")
pvals_wide <- experiments_to_wide(cors_df,values_col = "p.value") # Re-format data
pvals_wide_adj <- AdjustPmat(p_mat = as.matrix(pvals_wide[,names(exp_weights)])) # Adjust extreme pvalues 
path_pval <- CombinePval(pvals_wide_adj,exp_weights,n_cores = ncores)

ppval_df <- data.frame(Pathway.A=gsub(";.*","",names(path_pval)), # Re-format data 
                       Pathway.B=gsub(".*;","",names(path_pval)), 
                       p.value=path_pval) %>% 
            mutate(p.adj=p.adjust(p = p.value, method = "fdr")) # Adjust p-values for multiple testing
rownames(ppval_df)<-NULL
        
# Store partial table
write.table(ppval_df,file = file.path(res_folder,"pathway_pvalues.tsv"),
            sep="\t",col.names=TRUE,row.names=FALSE,quote=FALSE)

rm(pvals_wide,path_pval) # Clean up environment
      
# Join all tables 

pval_cor<-dplyr::inner_join(ppval_df,pcor_df,by=c("Pathway.A","Pathway.B"))
res.out<-dplyr::right_join(pval_cor,res,by=c("Pathway.A","Pathway.B"))
write.table(res.out,file = file.path(res_folder,"pathway_estimates.tsv"),
            sep="\t",col.names=TRUE,row.names=FALSE,quote=FALSE)

# Time elapsed
globalend <- Sys.time()
str<-globalend-globalstart
cat("\nFinished.\n")
print(str)

        