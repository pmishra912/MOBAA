

#' plot_feature_boxplot
#'
#' @description Plot boxplot to visualize distributional differences in feature expression/methylation/abundance between individuals in bicluster and those outside the bicluster.
#'
#' @param omics_list Named list of omic data matrices (features × samples).
#' @param result_list Output from analyze_multiomic_relation() (a list of data frames).
#' @param multiomic_relation multiomic_relation object
#' @param omic_type Name of the omic layer to plot (e.g., "expression", "methylation").
#' @param rank Which ranked feature to plot (e.g., 1 (default) for top feature by p-value).
#' @return Boxplot
#' @examples
#' plot_feature_boxplot()
#' @export
plot_feature_boxplot <- function(omics_list, result_list, multiomic_relation, omic_type, rank = 1) {
  
  
  if (!(omic_type %in% names(omics_list))) {
    stop("Specified omic_type not found in omics_list.")
  }
  if (!(omic_type %in% names(result_list))) {
    stop("Specified omic_type not found in result_list.")
  }

  omic_data <- omics_list[[omic_type]]
  result_df <- result_list[[omic_type]]

  if (nrow(result_df) < rank) {
    stop("Rank exceeds number of features available.")
  }

  # Get feature name at requested rank
  feature <- result_df$feature[rank]

  # Get expression values for the feature
  if (!(feature %in% rownames(omic_data))) {
    stop("Feature not found in omic data.")
  }

  values <- omic_data[feature, ]

  # Determine sample groups (clean names from multiomic_relation)
  cleaned_sample_names <- setNames(
    multiomic_relation$sample_ids,
    nm = sub("_.*$", "", names(multiomic_relation$sample_ids))
  )

  related_samples <- cleaned_sample_names[[omic_type]]
  all_samples <- colnames(omic_data)

  related_samples <- intersect(all_samples, related_samples)
  unrelated_samples <- setdiff(all_samples, related_samples)

  # Create data frame for plotting
  plot_df <- data.frame(
    sample = c(related_samples, unrelated_samples),
    value = as.numeric(values[c(related_samples, unrelated_samples)]),
    group = c(rep("related", length(related_samples)),
              rep("unrelated", length(unrelated_samples)))
  )

  # Basic boxplot using ggplot2
  library(ggplot2)

  p <- ggplot(plot_df, aes(x = group, y = value, fill = group)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.6) +
    geom_jitter(width = 0.2, alpha = 0.8) +
    labs(
      title = paste("Feature:", feature, "| Omic:", omic_type),
      y = "Expression / Omic Value",
      x = "Sample Group"
    ) +
    theme_minimal() +
    theme(legend.position = "none")
    print(p)
}


