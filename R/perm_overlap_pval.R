
#' perm.overlap.sig
#'
#' @description Calculates permutation based p-value for sample overlap among the biclusters from the different omics.
#'
#' @param omic_relations Output from 'find_relations'
#' @param sample_list Vector of sample names in the omics data (all omics must have exact same sample list).
#' @return P-value from prop.test
#' @examples
#' risk_score_relation()
#' @export

perm.overlap.sig <- function(omic_relations, sample_list){
    
    overlap_rand <- rep(NA,10000) # random overlap of samples of same size as the BC
    
    for(i in 1:10000){
        
        omic_rand_samples <- vector("list",length=length(omic_relations))
        
        for(j in 1:length(omic_relations)){
            omic_rand_samples[[j]] <- sample(sample_list,length(omic_relations[[j]]))
        }
        
        overlap_rand[i] <- length(Reduce(intersect, omic_rand_samples))
    }
    
    samples_overlap <- length(Reduce(intersect, omic_relations))
    pval <- (length(which(overlap_rand>=samples_overlap))+1)/10000
    return(pval)
}

