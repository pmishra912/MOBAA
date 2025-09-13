
#' correlate_eigenfeatures
#'
#' @description Computes a correlation matrix between the summary expression levels (PC1s) of different omic biclusters.
#'
#' @param omic_list A named list of omic data matrices (genes, CpGs, etc.). Each matrix has features (rows) and samples (columns).
#' @param multiomic_relation Multiomic bicluster relations. Output from 'find_relations'
#' @param output_prefix Prefix for output name. For example, "eigen_cor_plot"
#' @return List object with correlation matrix, p-value matrix and pairwise scatter plots.
#' @examples
#' correlate_eigenfeatures()
#' @export
correlate_eigenfeatures <- function(omic_list, multiomic_relation, output_prefix = "eigen_cor_plot") {
    
  if (!requireNamespace("ggpubr", quietly = TRUE)) stop("Package 'ggpubr' is required.")
  if (!requireNamespace("factoextra", quietly = TRUE)) stop("Package 'factoextra' is required.")

  library(ggpubr)
  library(factoextra)

  # STEP 1: Get overlapping samples across all omics
  shared_samples <- Reduce(intersect, multiomic_relation$sample_ids)
  message("Number of overlapping samples: ", length(shared_samples))

  if (length(shared_samples) < 3) {
    warning("Too few overlapping samples for PCA/correlation. Returning NULL.")
    return(NULL)
  }

  eigenfeatures <- list()

  for (omic in names(multiomic_relation$features)) {
    # Check if omic exists in omic_list
    if (!omic %in% names(omic_list)) {
      warning(paste("Omic", omic, "not found in omic_list. Skipping."))
      next
    }

    omic_matrix <- omic_list[[omic]]
    selected_features <- multiomic_relation$features[[omic]]

    # Subset matrix
    sub_matrix <- omic_matrix[
      rownames(omic_matrix) %in% selected_features,
      colnames(omic_matrix) %in% shared_samples,
      drop = FALSE
    ]

    if (nrow(sub_matrix) < 2) {
      warning(paste("Not enough features in omic", omic, "for PCA. Skipping."))
      next
    }

    # PCA (samples as rows)
    pca_result <- prcomp(t(sub_matrix), scale. = TRUE)
    pc1_scores <- pca_result$x[, 1]
    eigenfeatures[[omic]] <- pc1_scores
  }

  # Combine PC1s
  eigen_df <- as.data.frame(eigenfeatures)

  # Correlation matrix
  cor_matrix <- cor(eigen_df, use = "pairwise.complete.obs")

  # P-value matrix
  omic_names <- names(eigen_df)
  pval_matrix <- matrix(NA, nrow = length(omic_names), ncol = length(omic_names),
                        dimnames = list(omic_names, omic_names))

  for (i in seq_along(omic_names)) {
    for (j in seq_along(omic_names)) {
      if (i < j) {
        test <- cor.test(eigen_df[[i]], eigen_df[[j]], method = "pearson")
        pval_matrix[i, j] <- test$p.value
        pval_matrix[j, i] <- test$p.value
      } else if (i == j) {
        pval_matrix[i, j] <- 0
      }
    }
  }

  # Pairwise plots
  plot_list <- list()
  combos <- combn(names(eigen_df), 2, simplify = FALSE)

  for (pair in combos) {
    df_pair <- data.frame(
      PC1_1 = eigen_df[[pair[1]]],
      PC1_2 = eigen_df[[pair[2]]]
    )

    p <- ggscatter(
      df_pair, x = "PC1_1", y = "PC1_2",
      add = "reg.line", conf.int = TRUE,
      cor.coef = TRUE, cor.method = "pearson",
      xlab = paste(pair[1], "PC1"),
      ylab = paste(pair[2], "PC1"),
      #title = paste("Correlation between", pair[1], "and", pair[2])
    )

    filename <- paste0(output_prefix, "_", pair[1], "_vs_", pair[2], ".pdf")
    ggsave(filename, plot = p, width = 6, height = 5)

    plot_list[[paste(pair, collapse = "_vs_")]] <- p
  }

  return(list(
    correlation_matrix = cor_matrix,
    p_value_matrix = pval_matrix,
    plots = plot_list,
    eigenfeatures = eigen_df
  ))
}
