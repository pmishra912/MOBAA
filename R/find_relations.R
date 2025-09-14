
#' find_relations
#'
#' @description Finds all possible combinations across the input omic biclusters, calculates overlaps for each multi-omic set and return the best overlapping combination.
#'
#' @param biclusters_list A named list of biclusters for all omics. Output from 'find_ensemble_bc'. Must be named as c("gwe","meth","prot").
#' @param module Bicluster module name
#' @param pheno Data frame with sample IDs as one of the columns with name 'sampleID'
#' @param bc_dist_obj Output from 'generate_jaccard_dist' for bicluster matrix and jaccard dist for sample ids and bicluster names respectively
#' @param bc_hclust_obj Output from 'hclust_bcs' function for module names/colours
#' @return The best multi-omic relations based on sample overlap count.
#' @examples
#' find_relations()
#' @export
find_relations <- function(biclusters_list,module,pheno,bc_dist_obj,bc_hclust_obj) {
    
    ############# generate named list of samples from each biclusters in a module #########
    
    bc_dist_mat <- as.matrix(bc_dist_obj[[1]])
    bc_matrix <- bc_dist_obj[[2]]
    
    rownames(bc_matrix) <- pheno$sampleID # attach samples names to the bc matrix. The order should be checked beforehand.
    
    module_assignment <- bc_hclust_obj[[3]]
    module_assignment <- module_assignment[match(colnames(bc_dist_mat),module_assignment$biclusters),]
    
    if(sum(module_assignment$biclusters==colnames(bc_matrix))<ncol(bc_matrix)){
        print("BC names in hclust object and dist object do not match!")
    }
    
    
    module_bcs <- colnames(bc_dist_mat)[module_assignment$modules==module]

    # extract samples
    bc_samples <- list()


    for(i in 1:length(module_bcs)){
        ind <- which(colnames(bc_matrix)==module_bcs[i])
        bc_samples[[i]] <- rownames(bc_matrix)[which(bc_matrix[,ind]==TRUE)]
    }

    names(bc_samples) <- module_bcs
    
    
    ################################ find relations ##############################
    
    ## Get unique omics prefixes (assumes names are like "omic_bcid")
    omics <- unique(unlist(lapply(strsplit(names(bc_samples), "_"), function(x) x[[1]])))

    # Create a list of vectors for each omic
    vector.list <- lapply(omics, function(omic) {
        bc_samples[grep(paste0("^", omic, "_"), names(bc_samples))]
    })

    # Generate all combinations across all omics
    bc_names_list <- lapply(vector.list, names)
    all.possible.comb <- expand.grid(bc_names_list, stringsAsFactors = FALSE)

    # Calculate overlaps for each combination
    overlap_count <- apply(all.possible.comb, 1, function(row) {
        selected_bcs <- bc_samples[match(row, names(bc_samples))]
        length(Reduce(intersect, selected_bcs))
    })

    # Find the combination with the maximum overlap
    best_comb_index <- which.max(overlap_count)
    best_comb_names <- all.possible.comb[best_comb_index, ]

    # Extract corresponding entries from bc_samples
    selected_indices <- match(unlist(best_comb_names), names(bc_samples))
    omic_relation <- bc_samples[selected_indices]

    
    ###################### extract features from the identified relations #######################
    # 'omic_relation' is a named list of sample ids from multiomic BCs in a relation. The name is in form 'omic_type_number'.
    # Extract the numeric part of the name (e.g., 3 from "gwe_3"). Use that number to extract the appropriate bicluster from biclusters. Extract @rowname from the bicluster.
    
    # Extract numeric part after the underscore
    bic_index <- as.numeric(sub(".*_", "", names(omic_relation)))
    bic_omic_names <- sub("_.*", "", names(omic_relation))
    
    
    
    # Extract the biclusters for each omic
    features_list <- list()
    
    for(i in 1:length(names(omic_relation))){
        bc_tmp <- biclusters_list[[match(bic_omic_names[i],names(biclusters_list))]] # matches the bic_omic_names
        bic <- bc_tmp[[bic_index[i]]]
        
        features <- tryCatch({
          bic@rowname
        }, error = function(e) {
          warning(paste("Failed to extract @rowname from bicluster", names(omic_relation)[i], "-", e$message))
          return(NULL)
        })
        features_list[[i]] <- features
    }
    names(features_list) <- bic_omic_names
 
 relations <- list(sample_ids = omic_relation, features = features_list)
 return(relations)
 
}
