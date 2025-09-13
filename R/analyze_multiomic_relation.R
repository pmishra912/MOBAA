
#' analyze_multiomic_relation
#'
#' @description Performs a two-group test on each feature between bicluster and out-of-bicluster samples.
#'
#' @param omics_list A named list of omics matrices (each matrix: rows = features, columns = samples).
#' @param multiomic_relation multiomic_relation object with sample_ids: list mapping omic names to vector of sample IDs in biclusters and features: list mapping omic names to vector of features in biclusters.
#' @param test A test type: "t-test" (default) or "wilcox".
#' @param inter.omics.cor Logical. TRUE if feature-by-feature correlation is desired.
#' @param pval_threshold Default is 0.05. Only needed when inter.omics.cor is TRUE
#' @return A data.frame with feature name, omic name, estimate (mean difference), p-value.
#' @examples
#' analyze_multiomic_relation()
#' @export
analyze_multiomic_relation <- function(omics_list, multiomic_relation, test = "t-test", inter.omics.cor = FALSE, pval_threshold = 0.05) {
  
  # inter.omics.cor: logical, TRUE if feature-by-feature correlation is needed,
  # pval_threshold: default is 0.05. Only needed when inter.omics.cor is TRUE
  
  if (!(test %in% c("t-test", "wilcox"))) {
    stop("Invalid test specified. Choose 't-test' or 'wilcox'.")
  }

  # Validate input
  if (is.null(names(multiomic_relation$features))) {
    stop("multiomic_relation$features must be a named list.")
  }
  if (is.null(names(multiomic_relation$sample_ids))) {
    stop("multiomic_relation$sample_ids must be a named list.")
  }

  # Clean names
  cleaned_feature_names <- setNames(
    multiomic_relation$features,
    nm = sub("_.*$", "", names(multiomic_relation$features))
  )
  cleaned_sample_names <- setNames(
    multiomic_relation$sample_ids,
    nm = sub("_.*$", "", names(multiomic_relation$sample_ids))
  )

  final_results <- list()
  subsetted_omics <- list()

  # --- 1. Differential Test for Each Omic ---
  for (omic_name in names(omics_list)) {
    omic_data <- omics_list[[omic_name]]
    omic_results <- list()

    if (!(omic_name %in% names(cleaned_feature_names)) ||
        !(omic_name %in% names(cleaned_sample_names))) {
      next
    }

    related_samples <- cleaned_sample_names[[omic_name]]
    related_features <- cleaned_feature_names[[omic_name]]

    common_features <- intersect(rownames(omic_data), related_features)
    if (length(common_features) == 0) next
    omic_data <- omic_data[common_features, , drop = FALSE]

    all_samples <- colnames(omic_data)
    unrelated_samples <- setdiff(all_samples, related_samples)
    related_samples <- intersect(all_samples, related_samples)

    # Save subsetted data for correlation
    subsetted_omics[[omic_name]] <- omic_data[common_features, related_samples, drop = FALSE]

    if (length(related_samples) < 2 || length(unrelated_samples) < 2) {
      warning(paste("Too few samples in one of the groups for omic:", omic_name))
      next
    }

    for (feature in common_features) {
      values <- omic_data[feature, ]
      group1 <- as.numeric(values[related_samples])
      group2 <- as.numeric(values[unrelated_samples])

      test_estimate <- NA
      direction_estimate <- NA
      p_value <- NA

      if (test == "t-test") {
        test_result <- tryCatch(t.test(group1, group2), error = function(e) NULL)
        if (!is.null(test_result)) {
          test_estimate <- unname(diff(test_result$estimate))
          p_value <- test_result$p.value
        }
        direction_estimate <- mean(group1, na.rm = TRUE) - mean(group2, na.rm = TRUE)
      } else {
        test_result <- tryCatch(
          wilcox.test(group1, group2, exact = FALSE, conf.int = TRUE),
          error = function(e) NULL
        )
        if (!is.null(test_result) && !is.null(test_result$estimate)) {
          test_estimate <- test_result$estimate
          p_value <- test_result$p.value
        }
        direction_estimate <- mean(group1, na.rm = TRUE) - mean(group2, na.rm = TRUE)
      }

      direction <- if (is.na(direction_estimate)) {
        NA
      } else if (direction_estimate > 0) {
        "upregulated"
      } else if (direction_estimate < 0) {
        "downregulated"
      } else {
        "no change"
      }

      if (!is.na(p_value)) {
        omic_results[[length(omic_results) + 1]] <- data.frame(
          feature = feature,
          test_estimate = test_estimate,
          p_value = p_value,
          direction = direction,
          stringsAsFactors = FALSE
        )
      }
    }

    if (length(omic_results) > 0) {
      omic_df <- do.call(rbind, omic_results)
      omic_df <- omic_df[order(omic_df$p_value), ]
      rownames(omic_df) <- NULL
      final_results[[omic_name]] <- omic_df
    } else {
      message(paste("No valid results for omic:", omic_name))
    }
  }

  # --- 2. Generalized Inter-Omics Correlation ---
  
  if(inter.omics.cor){
  correlation_results <- list()
  omic_names <- names(subsetted_omics)

  if (length(omic_names) >= 2) {
    omic_combinations <- combn(omic_names, 2, simplify = FALSE)

    for (pair in omic_combinations) {
      omic1 <- pair[1]
      omic2 <- pair[2]
      mat1 <- subsetted_omics[[omic1]]
      mat2 <- subsetted_omics[[omic2]]

      # Common samples for this pair
      common_samples <- intersect(colnames(mat1), colnames(mat2))
      if (length(common_samples) < 3) next

      mat1 <- mat1[, common_samples, drop = FALSE]
      mat2 <- mat2[, common_samples, drop = FALSE]

      for (f1 in rownames(mat1)) {
        for (f2 in rownames(mat2)) {
          x <- as.numeric(mat1[f1, ])
          y <- as.numeric(mat2[f2, ])
          if (all(is.na(x)) || all(is.na(y))) next

          corr_test <- tryCatch(cor.test(x, y, method = "pearson"), error = function(e) NULL)
          if (!is.null(corr_test) && !is.na(corr_test$p.value) && corr_test$p.value <= pval_threshold) {
            correlation_results[[length(correlation_results) + 1]] <- data.frame(
              omic1 = omic1,
              feature1 = f1,
              omic2 = omic2,
              feature2 = f2,
              correlation = corr_test$estimate,
              p_value = corr_test$p.value,
              stringsAsFactors = FALSE
            )
          }
        }
      }
    }
  }

  correlation_df <- if (length(correlation_results) > 0) {
    do.call(rbind, correlation_results)
  } else {
    data.frame()
  }

  return(list(
    differential_tests = final_results,
    inter_omics_correlation = correlation_df
  ))
  
 }
  
  else{
      return(final_results)
  }
}
