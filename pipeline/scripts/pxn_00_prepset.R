############################################################
# Script: pdxn_00_prepset.R
#
# Author: Larisa M. Soto

# Description:
# This script processess an input gene set for the pdxn pipeline
#
# Arguments:
# - indir: directory containing input gene set objects and background dataset
# - outdir: directory to store outputs
# - dataset_name: name of the background reference dataset
# - geneset_name: name of the gene set 
# - tissue_select: name of the tissue
# - ncores: number of cores 
# - rels_char: character string for relationships to calculate
# - pcor_choice: choice of partial correlation
# - min_samples: min number of samples per tissue
#
# Usage:
# Rscript pdxn_00_prepset.R ~/projects/pdxn_2.0/data/gene_sets/genedex_2024/pathway_table_raw.csv gtextoil_gfilter genedex_gtextoil_gfilter ../input/gene_sets/genedex_gtextoil_gfilter ../input/gene_expression/gtextoil_gfilter ../input/gene_expression/genes_gtextoil_gfilter.csv 500 20 85 8
############################################################

start <- Sys.time() # get start time

suppressPackageStartupMessages({
    library(dplyr)
    library(plyr)
    library(data.table)
    library(furrr)
    require(combinat)
    require(future)
    require(purrr)
})
source('funcs/misc.R')
options(stringsAsFactors=FALSE)
cmd_args <- commandArgs(trailingOnly = T)

# User inputs - script args

csvfile <- cmd_args[1] # CSV file with standard gene set table
dsname <- cmd_args[2] # Dataset name
gsname <- cmd_args[3] # Geneset name
gs_dir <- cmd_args[4] # Geneset subfolder/ output folder
dset_dir <- cmd_args[5] # Dataset folder/ input folder with pre-processed background dataset
gene_file <- cmd_args[6] # Gene universe file (generated during pre-processing) 
max_genes <- as.numeric(cmd_args[7]) # Largest pathway size allowed
min_genes <- as.numeric(cmd_args[8]) # Smallest pathway size allowed 
max_jacq <- as.numeric(cmd_args[9])/100 # Max jacquard index between two pathways
ncores <- as.numeric(cmd_args[10]) # Number of cores to use in parallel
ag_dir <- cmd_args[11] # Optional, directory with references used to augment gene sets 
colln <- cmd_args[12] # Optional, name of collection used to augment the custom geneset
ntype <- cmd_args[13] # Optional, type of network to calculate between custom gene set and collection: unipartite, bipartite, uni-bipartite
augment<-FALSE # Flage to control augment mode (FALSE by default)

###############################################
# Main
###############################################

message("\nINPUT ARGS\n")
message("GENESET TABLE = ",csvfile)
message("GENE SET NAME = ",gsname)
message("DATASET NAME = ",dsname)
message("GENESET DIRECTORY = ",gs_dir)
message("DATASET DIRECTORY = ",dset_dir)
message("GENE UNIVERSE = ",gene_file)
message("MAX SET SIZE = ",max_genes)
message("MIN SET SIZE = ",min_genes)
message("MAX JI OF PAIRS = ",max_jacq)
message("NUMBER OF CORES = ",ncores)

if(exists("colln") && !is.na(colln) && colln!='' && exists("ntype") && !is.na(ntype) && ntype!=''){
    # If augmenting gene set 
    message("AUGMENT DIRECTORY = ",ag_dir)
    message("COLLECTION = ",colln)
    message("NETWORK TYPE = ",ntype) 
    augment<-TRUE # Flage to enable augment mode
}
message('\n')

# Start parallel multisession

plan(multisession, workers = ncores)

# Prepare dataset structure for PDxN 
####################################################

create_directory(gs_dir)
out_csvfile <- file.path(gs_dir,"pathway_table.csv") # gene set table 
out_rdsfile <- file.path(gs_dir,"pathway_list.RDS") # final gene set rds object 
out_stat_file <- file.path(gs_dir,"pathway_pairs_stats") # pathway pairs stats table

# Load data
####################################################

message("Reading input gene set and gene universe...")
gs_df <- fread(csvfile)
genes_ref <- read.table(gene_file,header=TRUE)
message('    Total number of pathways = ',length(unique(gs_df$set_name)))
message('    Total number of genes in universe = ',length(unique(genes_ref$entrez_id)))

# Process gene set
####################################################

gs_df <- gs_df %>%
             mutate(genes = as.character(genes),
                    node_type = 'core') 

## Augmenting core gene set with collection

if(augment){
    # Load collection
    message('Loading collection...')
    cl_df<-load_collection(coll_name = colln,aug_dir = ag_dir) %>%
            mutate(genes = as.character(genes),
                   node_type = 'augmented')
    message('    Total number of pathways in collection: ',length(unique(cl_df$set_name)))

    message('Augmenting gene set with collection...')    
    gs_df <- rbind(gs_df,cl_df)    
    message('    Final number of pathways: ',length(unique(gs_df$set_name)))    
}

## Filter genes in gene universe 
message("Filtering genes and genesets...")

gs_df <- gs_df %>%
         filter(genes %in% as.character(genes_ref$entrez_id)) 
message('    Number of set genes within universe = ',length(unique(gs_df$genes)))

## Filter pathway sizes
pwy_set <- gs_df %>%
            as.data.frame() %>%
            dplyr::group_by(set_name) %>%
            dplyr::summarize(ngenes=length(unique(genes))) %>%
            filter(ngenes>min_genes, # Keep pathways larger than min_genes
                   ngenes<max_genes) # Keep pathways smaller than max_genes

gs_df_filt <- gs_df %>%
              filter(set_name %in% pwy_set$set_name)
message('    Number of pathways within size bounds = ',length(unique(gs_df_filt$set_name)))

## Write gene set list object 
message("Creating PDxN geneset object...")
res <- gs_table_to_rds(gs_df = gs_df_filt,outfile = out_rdsfile)

# Process pathway pairs
####################################################

## Calculate overlap between pathway pairs on post processed gene set
if(augment){
    message('Building pathway pairs of ',ntype,' network...')
    stats<-get_pathway_pair_stats(gs_df = gs_df_filt,nwk_type = ntype)
}else{
    message('Building pathway pairs of unipartite network...')
    stats<-get_pathway_pair_stats(gs_df = gs_df_filt) 
}

## Filter pathway pairs 
pwy_pairs <- stats %>%
             filter(Jacquard.Index<max_jacq) # Discard pairs with JI > max_jacq
message('    Number of pathway pairs to compare = ',nrow(pwy_pairs),' out of ',nrow(stats))

## Write tables
message('Writing outputs...')
write.table(pwy_pairs,
            file = paste0(out_stat_file,".csv"), 
            quote=T, row.names = F, col.names =T, sep = ",")

write.table(gs_df_filt,
            file=out_csvfile,
            quote=F,sep=",",row.names = F)


end <- Sys.time() # get start time
# print elapsed time
print(end - start) # calculate difference
message('Finished.')


