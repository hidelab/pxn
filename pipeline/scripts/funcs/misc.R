######################################################
# Create directory 
######################################################

create_directory <- function(dir_path) {
  # Check if the directory exists, if not, create it
  if (!dir.exists(dir_path)) {
    dir.create(dir_path, recursive = TRUE)
  }
}

######################################################
# Get expression data - deprecated  
######################################################

getExprs <- function(x,ext=".collapse.RDS",...){
    
# Arguments 
# x - tissue name 
# ext - object file extension after tissue name

    tissue_fn <- gsub("[#,%:]",".",x)
    # get path to sample of a given tissue
    tissue_rds <- paste0(barcode_dir,tissue_fn,"/",tissue_fn,ext)
    # load normalized expression values
    tissue_exprs <- readRDS(tissue_rds)$datETcollapsed
    return(tissue_exprs)
}

######################################################
# Generate run report
######################################################

generate_report <- function(geneset_file, r_cores, output_folder, pcor_choice, rels, number_of_pathways) {

# Example usage:
# report_df <- generate_report(geneset_file, r_cores, output_folder, pcor_choice, rels, number_of_pathways)


# Updates to the original code:
# - The function arguments are used directly without creating new variables.
# - The `rels_mapping` list is used to map the `rels` value to the corresponding text variant.
# - The rest of the function remains the same, ensuring that the report table is generated correctly as a data frame.

  # Initiate reporter matrix
  report <- matrix(nrow = 0, ncol = 2)

  # Adding rows to reporter matrix
  report <- rbind(report, c("Genesets file", geneset_file))
  report <- rbind(report, c("Genesets file creation Time and Date", format(file.info(geneset_file)$ctime, "%a %b %d %X %Y")))
  report <- rbind(report, c("Cores", r_cores))
  report <- rbind(report, c("Output folder", output_folder))
  report <- rbind(report, c("Partial Correlation", pcor_choice))
  
  # Adding the relation row using the rels_mapping list
  
  relation_text <- map_rels_code(rels)
    
  report <- rbind(report, c("Relation", paste(rels, relation_text, sep = " > ")))
  report <- rbind(report, c("Run Date & Time", format(Sys.time(), "%a %b %d %X %Y")))
  report <- rbind(report, c("Pairs", number_of_pathways))

  # Convert matrix to data frame
  report_df <- as.data.frame(report, stringsAsFactors = FALSE)
  colnames(report_df) <- c("Description", "Value")
  
  return(report_df)
}

######################################################
# Filter tissue matrices
######################################################

filter_tissue_expr_list <- function(cts_list,min_counts,min_samples,gu_file){

    # Filter the list of gene expression matrices on the gene level to only keep
    # genes with at least `min_counts` counts in at least `min_samples` samples.
    # Will write the gene universe into gu_file. 
    # Arguments:
    #   - cts.list = list of tissue expression matrices 
    #   - min_counts = minimum number of gene counts
    #   - min_samples = minimum number of samples with min_counts
    # Returns a list with slots: 
    #   - filtered_cts = filtered list of tissue expression matrices (same length 
    #                     as input
    #   - gene_universe = data frame of entrez ids of genes passing the filters

    # Find gene universe
    gene_univ<-lapply(cts_list,function(cts,...){
                        g<-rowSums(cts>min_counts)>min_samples
                        names<-as.character(rownames(cts)[g])
                        return(names)
                    }) %>%
                do.call(c,.) %>%
                sort() %>%
                unique()

    # Filter tissue matrices 
    cts_list_filt <- lapply(cts_list,function(cts,...){
        
                                cts_filt<-cts[rownames(cts)%in%gene_univ,]
                                missing_genes<-gene_univ[!gene_univ%in%rownames(cts)]

                                # Add zeros for genes in the universe that are not expressed in a given sample
                                if(length(missing_genes!=0)){
                                    cts_filt[missing_genes,]<-0
                                }     

                                cts_filt<-cts_filt[gene_univ,]
        
                                return(cts_filt)
                            })

    # Store file with gene universe 
    #message('Gene universe contains ',length(gene_univ),' genes')
    #message('Wrote gene universe to file ',gu_file)

    gene_univ<-data.frame(entrez_id=as.character(gene_univ))
    write.table(gene_univ,file=gu_file,sep='\t',row.names=F,quote=T,col.names=T)

    # Outputs
    #message('Returning filtered matrices and gene universe')
    outs<-list('filtered_cts' = cts_list_filt,
               'gene_universe' = gene_univ)
    
    return(outs) 
}

######################################################
# Generate input file structure 
######################################################

tissue_list_to_dirs<-function(cts_list,output_dir,meta_file){

    # This function receives a named list with gene expression matrices and
    # creates a subdirectory for every named item where it writes the matrix
    # as an RDS file and as a CSV table using the list item name as file name. 
    #
    # Arguments:
    # cts_list = named list of expression matrices 
    # output_dir = path to generate subdirectory structure
    # 
    # Returns: a message with how many directories where created and the size 
    # of each matrix
        
    out<-names(cts_list) %>%
         lapply(.,function(tissue,...){
        
            objdir<-file.path(output_dir,"subgroups",tissue)
            create_directory(objdir)
    
            cts<-cts_list[[tissue]]
            write.table(cts,file=file.path(objdir,paste0(tissue,".csv")),quote=TRUE,row.names=TRUE,col.names=TRUE,sep=",")
             
            return(dim(cts)) 
         })
    
    names(out)<-names(cts_list)

    # Create dataset metadata file
    
    meta<-data.frame(subgroup=names(out),
                     nsamples=lapply(out,function(l){l[[2]]}) %>% unlist())
    write.table(meta,file=meta_file,row.names = FALSE, quote=TRUE, sep=",")
    
    res<-list("ndirs"= length(out),
              "sizes"= out)

    return(res)
    
}


######################################################
# Generate input gene set object 
######################################################

gs_table_to_rds<-function(infile=NULL,gs_df=NULL,outfile){

    # This function receives a gene set table either as a file or as 
    # a data frame object and creates an RDS object with a named
    # gene set list. 
    # 
    # Arguments:
    # infile = input file with the gene set table with two columns the first
    # one is has the names of the gene sets and the second one the genes. 
    # gs_df = a data frame object with two columns, the first one has the gene
    # set names and the second one the genes.
    # outfile = file to write the RDS object with the gene set list
    # Returns:
    # A list with the function call and the number of sets in the object

    if( (is.null(gs_df) & is.null(infile)) | (!is.null(gs_df) & !is.null(infile)) ){
        stop("Error: must provide either a data frame or an input file")
    }
    
    if(is.null(gs_df)){
        gs_df<-fread(infile)
    }    
    colnames(gs_df)[1:2]<-c("set_name","genes")   
    gs_df<-gs_df %>%
            group_by(set_name) %>%
            dplyr::summarize(geneset=list(genes))

    gs_list<-gs_df$geneset
    names(gs_list)<-gs_df$set_name
    saveRDS(gs_list,outfile)

    res<-list("args"=list(infile,gs_df,outfile),
              "nsets"=length(gs_list))

    return(res)
}

######################################################
# Calculate pathway overlaps
######################################################

#' Calculate Pairwise Pathway Statistics
#'
#' This function computes pairwise statistics (e.g., overlap coefficient and Jaccard index)
#' for a set of pathways based on gene overlap. It supports unipartite and bipartite network configurations.
#'
#' @param gs_df A data frame containing pathway information with columns `set_name`, `node_type`, 
#'              and `genes`. The `genes` column should be a list or vector of gene identifiers for each pathway.
#' @param nwk_type A string specifying the type of network to create. Options are:
#'                 - `'unipartite'`: Compute statistics for all pairwise combinations of pathways.
#'                 - `'bipartite'`: Compute statistics for pathway pairs where one node is `core` 
#'                   and the other is `augmented`.
#'                 - `'uni_bipartite'`: Combine unipartite core-core pairs with bipartite core-augmented pairs.
#' @param ... Additional arguments to be passed to the `pwy_stats` function, if needed.
#'
#' @return A data frame with pairwise pathway statistics, including:
#'         - `Pathway.A`: The first pathway in the pair.
#'         - `Pathway.B`: The second pathway in the pair.
#'         - `Overlap.Coefficient`: The overlap coefficient between the two pathways.
#'         - `Jacquard.Index`: The Jaccard index between the two pathways.
#'         - `Size.A`: The number of genes in `Pathway.A`.
#'         - `Size.B`: The number of genes in `Pathway.B`.
#'
#' @import dplyr
#' @import furrr
#' @examples
#' # Example pathway data frame
#' gs_df <- data.frame(
#'   set_name = c("Pathway1", "Pathway2", "Pathway3"),
#'   node_type = c("core", "augmented", "core"),
#'   genes = list(c("Gene1", "Gene2"), c("Gene3", "Gene4"), c("Gene2", "Gene5"))
#' )
#'
#' # Calculate pairwise pathway statistics for a unipartite network
#' result <- get_pathway_pair_stats_2(gs_df, nwk_type = 'unipartite')
#'
#' @export

get_pathway_pair_stats <- function(gs_df,nwk_type='unipartite',...) {

  gs_df_list<-gs_df %>%
                group_by(set_name,node_type) %>%
                dplyr::summarize(geneset=list(genes),.groups='keep')
  ##message('    Reading ', nrow(gs_df_list), ' input pathways...')

  pathways <- gs_df_list$geneset
  names(pathways) <- gs_df_list$set_name
    
  ##message('    Creating table of pathway pairs...')
    
  if(nwk_type=='unipartite'){  
    # Everything vs everything
    ##message('    Building unipartite network with ',length(pathways),' nodes')
    pathway_pairs <- expand.grid(A = names(pathways), B = names(pathways)) %>%
                     filter(as.character(A) < as.character(B)) 
  }else{
      
    # Split nodes
    clist <- gs_df_list %>% filter(node_type=='core')
    core_nodes <- clist$set_name

    alist <- gs_df_list %>% filter(node_type=='augmented')
    aug_nodes <- alist$set_name
    
    # Bipartite network

    #message('    Building bipartite network...')
    pathway_pairs <- expand.grid(A = core_nodes, B = aug_nodes) 
    #message('        Number of pathways = ',length(core_nodes)+length(aug_nodes))
    #message('        Number of pairs = ',nrow(pathway_pairs))
        
    if(nwk_type=='uni_bipartite'){
        
        # Unipartite core network + bipartite
        #message('    Building uni_bipartite network...')
        uni_pairs <- expand.grid(A = core_nodes, B = core_nodes) %>%
                       filter(as.character(A) < as.character(B)) 
        pathway_pairs <- rbind(uni_pairs,pathway_pairs)
        #message('        Number of pathways = ',length(core_nodes)+length(aug_nodes))
        #message('        Number of pairs = ',nrow(pathway_pairs))
    }
      
    colnames(pathway_pairs) <- c('Pathway.A','Pathway.B')
  } 
    
  # Calculate pair stats
  #message('    Calculating pair statistics...')

  pwy_stats <- function(pathway_a, pathway_b,...) {
      
        genes_a <- unique(pathways[[pathway_a]])
        genes_b <- unique(pathways[[pathway_b]])
        
        intersection <- length(intersect(genes_a, genes_b))
        union <- length(genes_a) + length(genes_b) - intersection
        size_smallest <- min(length(genes_a), length(genes_b))
        
        overlap_coeff <- intersection / size_smallest
        jacq_index <- intersection / union
      
        outs<-data.frame(Pathway.A = pathway_a, 
                   Pathway.B = pathway_b,
                   Overlap.Coefficient = overlap_coeff,
                   Jacquard.Index = jacq_index,
                   Size.A = length(genes_a),
                   Size.B = length(genes_b))

      return(outs)
  }
    
  pathway_pairs_stats <- pathway_pairs   %>%
                         furrr::future_pmap_dfr(~ pwy_stats(..1, ..2)) # Apply pwy_stats for each pair
    
  return(pathway_pairs_stats)
}

######################################################
# Re-format all experiment-level into wide data frame
######################################################

#' Convert Experiments Data Frame to Wide Format
#'
#' This function converts a data frame containing experiment-level results from from long to wide format,
#' where the values from the specified column (estimates or p.values) are spread across different subgroups.
#'
#' @param df A data frame containing the experimental data. It must include the columns: "Pathway.A", "Pathway.B", "experiment", and the column specified by `values_col`.
#' @param values_col A character string specifying the name of the column containing the values to be spread across experiments. Default is "estimate".
#'
#' @return A data frame in wide format with the values spread across different experiments. The rows are named using a combination of "Pathway.A" and "Pathway.B".
#' @examples
#' library(dplyr)
#' library(tidyr)
#' library(tibble)
#'
#' # Sample data frame
#' df <- tibble(
#'   Pathway.A = c("A1", "A1", "A2", "A2"),
#'   Pathway.B = c("B1", "B2", "B1", "B2"),
#'   experiment = c("exp1", "exp1", "exp2", "exp2"),
#'   estimate = c(0.5, 0.6, 0.7, 0.8)
#' )
#'
#' # Convert to wide format
#' df_wide <- experiments_to_wide(df)
#' print(df_wide)
#'
#' @export
experiments_to_wide<-function(df,values_col = "estimate"){

    cols<-c("Pathway.A","Pathway.B","experiment",values_col)
    
    df.out<-df[,..cols] 
    colnames(df.out)[4]<-"values_column"
    
    df.out <- df.out %>%
                 tidyr::pivot_wider(id_cols = c(Pathway.A,Pathway.B),
                                    names_from = "experiment",
                                    values_from = "values_column") %>%
                 mutate(pair=paste(Pathway.A,Pathway.B,sep=";")) %>%
                 tibble::column_to_rownames("pair") %>%
                 select(-Pathway.A,-Pathway.B)
    
    return(df.out)
}



######################################################
# Load and Combine Gene Set Collections
######################################################

#' Load and Combine Gene Set Collections
#'
#' This function loads a collection of gene sets from a specified directory, checks for
#' file validity, and combines the data into a single data frame.
#'
#' @param coll_name A string specifying the name of the collection to load. This is used to 
#'        match files in the `collections` directory.
#' @param aug_dir A string specifying the base directory where the `collections` and 
#'        `tables` subdirectories are located.
#'
#' @return A data frame combining all the gene set tables from the collection.
#'         If some files in the collection are missing, they are skipped with a warning.
#' @examples
#' # Assuming `aug_dir` contains a valid collections directory:
#' # load_collection("example_collection", "/path/to/aug_dir")
#'
#' @import data.table
#' @import purrr
#' @export
load_collection <- function(coll_name,aug_dir){

  colln_path <- paste0(aug_dir,'/collections/',coll_name,'.txt')
    
  # Check if the file containing the list of file names exists
  if (!file.exists(colln_path)) {
    stop("The collection file does not exist.")
  }
  #message('    Found collection: ',coll_name)
  
  # Read the file list into a character vector
  file_names <- readLines(colln_path,warn = FALSE)
  file_names <- paste0(aug_dir,'/tables/',file_names,'_pathway_table.csv')
  #message('    Loading collection of ',length(file_names),' gene sets') 
    
  # Use map_dfr to read and combine all files into a single data frame
  combined_df <- purrr::map_dfr(file_names, ~ {
    if (file.exists(.x)) {
      data.table::fread(.x) # Use fread for fast reading
    } else {
      warning("The file '", .x, "' does not exist and will be skipped.")
      NULL # Return NULL for non-existing files
    }
  })
  
  return(combined_df)
}