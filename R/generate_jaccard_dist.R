

#' generate_jaccard_dist
#'
#' @description Calculates Jaccard distance between the biclusters across all omics.
#'
#' @param mosbi_outputs List of mosbi biclustering outputs for individual omics
#' @param omic_names Names of the omics in the same order (eg, c(meth,expr,protein))
#' @param omic_data List of omic data in same order.
#' @return Jaccard distance between the biclusters across all omics.
#' @examples
#' generate_jaccard_dist()
#' @export
generate_jaccard_dist <- function(mosbi_outputs,omic_names, omic_data,quiet = TRUE){
   
    if (!requireNamespace("proxy", quietly = TRUE)) {
        stop("The 'proxy' package is required but not installed.")
      }
    
    # Define a safe wrapper
      run_safe <- function(expr) {
        if (quiet) {
          suppressMessages(suppressWarnings(capture.output(res <- try(expr, silent = TRUE))))
        } else {
          res <- try(expr, silent = FALSE)
        }
        return(res)
      }
      
    bc_matrix_list <- list()
    for(j in 1:length(mosbi_outputs)){
        
        ncol <- length(mosbi_outputs[[j]]) # num of cols for each BC
        nrow <- ncol(omic_data[[j]]) # col is samples names in omic data
        bc_mat <- matrix(NA,nrow,ncol)
        rownames(bc_mat) <- colnames(omic_data[[j]])
        colnames(bc_mat) <- paste(omic_names[j],c(1:ncol),sep="_")
        
        # the following makes 0/1 matrix for each BC/column
        for(i in 1:ncol){
            bc_mat[,i] <- rownames(bc_mat) %in% mosbi_outputs[[j]][[i]]@colname
        }
        bc_matrix_list[[j]] <- bc_mat
    }
    
    # combine all sample x BC matrices from different omics
    combined_bc_mat <- do.call("cbind",bc_matrix_list)
    bc_dist <- run_safe(proxy::dist(combined_bc_mat, by_rows = FALSE, method = "Jaccard"))
    res <- list(bc_dist=bc_dist,bc_matrix=combined_bc_mat)
    return(res)
}
