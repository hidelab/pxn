######################################################
# Compute overlap coefficients
######################################################

OverlapCoefficient <- function(x,y){
    # function to calculate the overlap coefficient between x and y
    # which is defined as the size of the intersection divided by the
    # size of the smaller set
    #
    # Args
    #   x: a vector
    #   y: a vector
    #
    # Returns
    #   the overlap coefficient, a number between 0 and 1
    
    length(intersect(x,y))/min(length(unique(x)),length(unique(y)))
}

######################################################
# Pathway level summary stats
######################################################

GetSummary <- function(dat,gs,sum_fun){
    # function to calculate the summary statistic for the pathway
    #
    # Args.
    #   dat: genes by samples matrix
    #   gs: vector with the names of the genes in the gene set
    #   sum_fun: function to calculate the summary
    #
    # Returns
    #   a 1 by samples vector with the summary statistic for the pathway
    
    if(length(gs) > 1){
        # calculate summary for pathways with more than 1 element
        return(sum_fun(dat[rownames(dat) %in% gs,]))
    }else{
        # return actual value for pathways with a single element
        return(dat[rownames(dat) %in% gs,])
    }
}

######################################################
# Correlation wrapper
######################################################

ShrinkCor <- function(x,y,method="pearson"){
    # wrapper to estimate the correlation coefficient between x and y using the 
    # shrinkage estimator from Schaffer & Strimmer (2005) [corpcor package]
    # and the corresponding t-statistic and p-value
    #
    # Args
    #   x: a vector with n observations
    #   y: a vector with n observations
    #   method: character to pick either the Pearson or Spearman correlation coefficient
    #
    # Returns
    #   a named vector with the correlation estimate, the sample size n, the t-statistic
    #   and its corresponding p-value
    
    # function to get the t-statistic
    GetStatistic <- function(r,n){r*sqrt((n-2)/(1-r^2))}
    # get sample size
    if(length(x) == length(y)){
        n <- length(x)
    }else{
        #cat("\n x and y have different lengths! >=( \n")
        return(NA)
    }
    # determine method
    selected_method <- match(method,c("pearson","spearman"))
    # Pearson correlation
    if(selected_method == 1){
        estimate <- cor.shrink(cbind(x,y),verbose=F)
        statistic <- GetStatistic(estimate[2,1],n)
        p.value <- 2*pt(-abs(statistic),n-2)
    }else if(selected_method == 2){
        estimate <- cor.shrink(cbind(rank(x),rank(y)),verbose=F)
        statistic <- GetStatistic(estimate[2,1],n)
        p.value <- 2*pt(-abs(statistic),n-2)
    }else{
        #cat("invalid method! >=( \n")
        return(NA)
    }
    # prepare results
    res <- c(estimate[2,1],n,statistic,p.value)
    names(res) <- c("estimate","n","statistic","p.value")
    return(res)
}

######################################################
# Partial correlation wrapper
######################################################

ShrinkPCor <- function(x,y,z,method="pearson"){
    # wrapper to estimate the partial correlation coefficient x,y|z using the 
    # shrinkage estimator from Schaffer & Strimmer (2005) [corpcor package]
    # and the corresponding t-statistic and p-value
    #
    # Args
    #   x: a vector with n observations
    #   y: a vector with n observations
    #   z: a vector with n observations
    #   method: character to pick either the Pearson or Spearman partial correlation coefficient
    #
    # Returns
    #   a named vector with the partial correlation estimate, the sample size n, the t-statistic
    #   and its corresponding p-value
    
    # function to get the t-statistic
    GetStatistic <- function(r,n){r*sqrt((n-3)/(1-r^2))}
    # get sample size
    if(length(x) == length(y) & length(z) == length(x)){
        n <- length(x)
    }else{
        
        #cat("x,y and z have different lengths! >=( \n")
        return(NA)
    }
    # determine method
    selected_method <- match(method,c("pearson","spearman"))
    # Pearson correlation
    if(selected_method == 1){
        cor.xyz <- cor.shrink(cbind(x,y,z),verbose=F)
        estimate <- cor2pcor(cor.xyz)[1,2] 
        statistic <- GetStatistic(estimate,n)
        p.value <- 2*pt(-abs(statistic),n-3)
    }else if(selected_method == 2){
        cor.xyz <- cor.shrink(cbind(rank(x),rank(y),rank(z)),verbose=F)
        estimate <- cor2pcor(cor.xyz)[1,2] 
        statistic <- GetStatistic(estimate,n)
        p.value <- 2*pt(-abs(statistic),n-3)
    }else{
        #cat("invalid method! >=( \n")
        return(NA)
    }
    # prepare results
    res <- c(estimate,n,statistic,p.value)
    names(res) <- c("estimate","n","statistic","p.value")
    return(res)
}

######################################################
# Check for specific interactions
######################################################

check_rel <- function(n1,n2,rel){
    switch(rel,
           "1"  = if((startsWith(n1, "Pathway.") & startsWith(n2, "CMAP.")) | (startsWith(n2, "Pathway.") & startsWith(n1, "CMAP."))){return(TRUE)}else {return(FALSE)},
           "2"  = if((startsWith(n1, "Pathway.") & startsWith(n2, "CTD.")) | (startsWith(n2, "Pathway.") & startsWith(n1, "CTD."))){return(TRUE)}else {return(FALSE)},
           "3"  = if((startsWith(n1, "Pathway.") & startsWith(n2, "PharmGKB.")) | (startsWith(n2, "Pathway.") & startsWith(n1, "PharmGKB."))){return(TRUE)}else {return(FALSE)},
           "4"  = if((startsWith(n1, "CMAP.") & startsWith(n2, "CTD.")) | (startsWith(n2, "CMAP.") & startsWith(n1, "CTD."))){return(TRUE)}else {return(FALSE)},
           "5"  = if((startsWith(n1, "CMAP.") & startsWith(n2, "PharmGKB.")) | (startsWith(n2, "CMAP.") & startsWith(n1, "PharmGKB."))){return(TRUE)}else {return(FALSE)},
           "6"  = if((startsWith(n1, "CMAP.up.") & startsWith(n2, "CMAP.down")) | (startsWith(n2, "CMAP.up") & startsWith(n1, "CMAP.down"))){return(TRUE)}else {return(FALSE)},
           "7"  = if((startsWith(n1, "Pathway.") & startsWith(n2, "L1000CDS2.")) | (startsWith(n2, "Pathway.") & startsWith(n1, "L1000CDS2."))){return(TRUE)}else {return(FALSE)},
           "8"  = if((startsWith(n1, "L1000CDS2.up.") & startsWith(n2, "L1000CDS2.down")) | (startsWith(n2, "L1000CDS2.up") & startsWith(n1, "L1000CDS2.down"))){return(TRUE)}else {return(FALSE)},
           "9"  = if((startsWith(n1, "L1000CDS2.up") & startsWith(n2, "L1000CDS2.up"))){return(TRUE)}else {return(FALSE)},
           "10" = if((startsWith(n1, "L1000CDS2.down") & startsWith(n2, "L1000CDS2.down"))){return(TRUE)}else {return(FALSE)},
	       "11" = if((startsWith(n1, "Pathway.") & startsWith(n2, "Pathway."))){return(TRUE)}else {return(FALSE)},
           "12" = if((startsWith(n1, "mirSets.") & startsWith(n2, "Pathway.")) | (startsWith(n1, "Pathway.") & startsWith(n2, "mirSets."))){return(TRUE)}else {return(FALSE)}
    )
}

######################################################
# Experiment-level estimates of a gene set pair
######################################################

ProcessPair<- function(pathway_a,pathway_b,...){
 
    # Retreive pathway gene sets
    gsA <- as.character(gs_lst[[pathway_a]])
    gsB <- as.character(gs_lst[[pathway_b]])
    
    # Retrieve pathway summaries 
    psA <- summary_dis_list[[pathway_a]]
    psB <- summary_dis_list[[pathway_b]]

    # shared genes
    gsAB <- intersect(gsA,gsB) 

    
    # Get correlation between the summaries for the unique genes
    if(pcor_choice == "0" & length(gsAB) > 0) {
        # If pathways share genes, estimate conditional correlation (on shared genes)
        summaryAB <- GetSummary(dat=exprs_rnk, gs=gsAB, colMeans)
        outs <- ShrinkPCor(x=psA,y=psB,z=summaryAB,method = "pearson")
    }else{
        # Otherwise estimate correlation between gene sets
        outs <- ShrinkCor(x=psA,y=psB,method = "pearson")
    }

    outs <- c('Pathway.A'=pathway_a,
              'Pathway.B'=pathway_b,
              outs)
    
    return(outs)
    
}

######################################################
# Experiment-level estimates of a gene set pair
######################################################

ProcessElement <- function(ic,...){

    # helper function to get the experiment-level estimates for a gene-set pair
    # The function below is intended to transform a 1-d list into a 2-d index list.
    # Note that the list is upper triangular. So if you want to calculate a specific
    # type of relationship (like the check rell above) you have to inspect 
    # in correct order
    
    i = ceiling((sqrt(8*(ic+1)-7)+1)/2)
    j = ic-choose(floor(1/2+sqrt(2*ic)),2)
    
    # pathway gene sets
    gsA=gs_lst[[i]]
    gsB=gs_lst[[j]]
    
    n1 <- names(gs_lst[i])
    n2 <- names(gs_lst[j])
    
    pass <- FALSE
     
    for (r in 1:length(rels)) {
        if(check_rel(n1,n2,rels[r])) {
            pass <- TRUE
        }
    }
    
    # Check if this pairs passes all relationship checks
    if(rels[1] == 666) pass <- TRUE
    if(!pass) return(NULL)

    # shared genes
    gsAB <- intersect(gsA,gsB)
    
    # get correlation between the summaries for the unique genes
    tmp = data.frame(Pathway.A=names(gs_lst)[i],Pathway.B=names(gs_lst)[j])
    
    if(pcor_choice == "0") {
        if(length(gsAB) > 0){
            # if pathways share genes, estimate conditional correlation (on shared genes)
            summaryAB = GetSummary(dat=exprs_rnk,gs=gsAB,colMeans)
            
            tmp = c(tmp,ShrinkPCor(
                x=unlist(summary_dis_list[[i]]),
                y=unlist(summary_dis_list[[j]]),
                z=summaryAB,
                method = "pearson"
            ))
        }else{
            # otherwise, estimate correlation between gene sets
            tmp = c(tmp,ShrinkCor(
                x=unlist(summary_dis_list[[i]]),
                y=unlist(summary_dis_list[[j]]),
                method = "pearson"
            ))
        }
    }else {
        tmp = c(tmp,ShrinkCor(
            x=unlist(summary_dis_list[[i]]),
            y=unlist(summary_dis_list[[j]]),
            method = "pearson"
        ))
    }
    
    # calculate overlap coefficient
    tmp$Overlap.Coeff= OverlapCoefficient(as.numeric(gs_lst[[i]]),as.numeric(gs_lst[[j]]))
    
    setTxtProgressBar(pb,ic)
    return(tmp)
    
}

######################################################
# Check pathway - pathway interaction
######################################################

check_path_path <- function(x){

    # Checks if we have a pathway - pathway interaction: TRUE or FALSE

    str1 <- unlist(strsplit(unlist(x), split = " - "))[1]
    str2 <- unlist(strsplit(unlist(x), split = " - "))[2]
    
    if(startsWith(str1, "Pathway.") & startsWith(str2, "Pathway.")){
        return(TRUE)
    }else {
        return(FALSE)
    }
}

######################################################
# Adjust p-values and correlation estimates
######################################################

AdjustPmat <- function(p_mat,eps=1e-16){

# adjust p-values and correlation estimates,
# otherwise the functions to combine the p-values
# cannot handle values near 0 and values near 1

    #res <- t(apply(p_mat,1,function(pval){
    #                            pval[pval <= eps] = eps
    #                            pval[pval >= 1-eps] = 1 - eps
    #                            return(pval)
    #                        }
    #              )
    #        )
    #return(res)

    # Modified by Larisa M. Soto 

    p_mat[p_mat<=eps]<-eps
    p_mat[p_mat>= (1-eps)]<-1-eps
    
    return(p_mat)
}

######################################################
# Get correlation estimates
######################################################

GetCorEstimate <- function(ic){
    
    # load experiment-level estimates
    my_rds = paste0(output_folder,"/mean_pcor2_barcode/",fname[ic],"_cpad_pathcor.RDS")
    myLst = readRDS(my_rds)[ pairs_chunks[[ cp ]] ]
    # extract correlation estimates
    tmp = unlist( sapply(myLst, function(x){x[[ "estimate" ]]}), use.names = F )
    setTxtProgressBar(pb,ic)
    
    return( tmp )
    
}

######################################################
# Get P-values
######################################################

GetPvals <- function(ic){
    
    # load experiment-level estimates
    my_rds = paste0(output_folder,"/mean_pcor2_barcode/",fname[ic],"_cpad_pathcor.RDS")
    myLst = readRDS(my_rds)[ pairs_chunks[[ cp ]] ]
    # extract p-value
    tmp = unlist( sapply(myLst, function(x){x[[ "p.value" ]]}), use.names = F )
    setTxtProgressBar(pb,ic)
    
    return( tmp )
}

######################################################
# Combine P-values
######################################################

#' Combine P-values in Parallel
#'
#' This function combines p-values from a matrix using weighted Z-tests, performed in parallel to speed up computation.
#'
#' @param pvals_mat A numeric matrix where each row represents a set of p-values to be combined.
#' @param w_vec A numeric vector of weights corresponding to the p-values.
#' @param n_cores An integer specifying the number of CPU cores to use for parallel processing. Defaults to the number of available cores.
#'
#' @return A named numeric vector of combined p-values, with names corresponding to the row names of `pvals_mat`.
#'
#' @details The function uses the `sumz` function from the `metap` package to combine p-values using a weighted Z-test. The computation is parallelized using the `mclapply` function from the `parallel` package, which can significantly reduce computation time for large matrices.
#'
#' @examples
#' \dontrun{
#' library(metap)
#' 
#' # Example p-values matrix
#' pvals_mat <- matrix(runif(100), nrow = 10)
#' 
#' # Example weights vector
#' w_vec <- runif(10)
#' 
#' # Combine p-values in parallel using all available cores
#' result <- CombinePval(pvals_mat, w_vec)
#' 
#' # Print the result
#' print(result)
#' }
#'
#' @export
CombinePval <- function(pvals_mat, w_vec, n_cores = detectCores()) {
    
    # Define a function to process each row
    process_row <- function(prow) {
        if (any(is.na(prow))) {
            return(NA)
        } else {
            return(sumz(p = prow, weights = w_vec)$p)
        }
    }
    
    # Use mclapply for parallel processing
    cp <- mclapply(1:nrow(pvals_mat), function(i) {
        process_row(pvals_mat[i, ])
    }, mc.cores = n_cores)
    
    cp <- unlist(cp)
    names(cp) <- rownames(pvals_mat)
    
    return(cp)
}
