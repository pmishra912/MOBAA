#' hclust_bcs
#'
#' @description Performs hierarchical clustering of multi-omic biclusters using Jaccard index based distance
#'
#' @param jaccard_dist This is analogous to dissTOM in standard WGCNA pipeline where columns are biclusters
#' @param plots Logical (default: TRUE) — whether to generate dendrogram plot
#' @param cutHeight Cut the tree at the height of distance (default: 0.8)
#' @param minModuleSize Minimum number of biclusters in a module (default: 3)
#' @param quiet Logical (default: TRUE) — whether to suppress warnings and messages
#'
#' @return A list: all modules, only multiomic modules, module assignments
#' @examples
#' hclust_bcs(jaccard_dist)
#' @export
hclust_bcs <- function(jaccard_dist, plots = TRUE, cutHeight = 0.8, minModuleSize = 3, quiet = TRUE) {

  if (!requireNamespace("WGCNA", quietly = TRUE)) {
    stop("The 'WGCNA' package is required but not installed.")
  }

  if (!requireNamespace("dynamicTreeCut", quietly = TRUE)) {
    stop("The 'dynamicTreeCut' package is required but not installed.")
  }
  
  # Optional message/warning suppression
  quiet_wrap <- function(expr) {
    if (quiet) {
      suppressMessages(suppressWarnings(expr))
    } else {
      expr
    }
  }

  # Compute dendrogram and distance matrix
  bcTree <- quiet_wrap(hclust(jaccard_dist, method = "average"))
  bc_dist_mat <- quiet_wrap(as.matrix(jaccard_dist))

  # Dynamic tree cut
  dynamicMods <- quiet_wrap(dynamicTreeCut::cutreeDynamic(
    dendro = bcTree,
    distM = bc_dist_mat,
    deepSplit = 2,
    pamRespectsDendro = FALSE,
    minClusterSize = minModuleSize,
    cutHeight = cutHeight
  ))

  # Convert module numbers to colors
  dynamicColors <- quiet_wrap(WGCNA::labels2colors(dynamicMods))

  # List of module names (colors)
  module_names <- names(table(dynamicColors))

  # Map each module to its member biclusters
  module_bcs <- lapply(module_names, function(mod) {
    colnames(bc_dist_mat)[dynamicColors == mod]
  })
  names(module_bcs) <- module_names

  # Identify multiomic modules
  bc_multiomic_yesno <- sapply(module_bcs, function(x) {
    omic_tags <- sapply(strsplit(x, "_"), function(y) y[[1]])
    length(unique(omic_tags)) > 1
  })

  multiomic_module_names <- module_names[bc_multiomic_yesno]
  module_bcs_multiomics <- module_bcs[multiomic_module_names]

  message("Modules with multiomic biclusters:\n", paste(multiomic_module_names, collapse = "\n"))

  # Module assignment per bicluster
  bc_module_assignment <- data.frame(
    biclusters = colnames(bc_dist_mat),
    modules = dynamicColors,
    stringsAsFactors = FALSE
  )

  # Plot dendrogram + colors
  if (plots && !is.null(bcTree) && !is.null(dynamicColors)) {
    if (length(dynamicColors) == length(bcTree$order)) {
      quiet_wrap({
          # Save to PDF
        pdf("bicluster_dendogram.pdf")
        WGCNA::plotDendroAndColors(
          bcTree,
          dynamicColors,
          "Dynamic Tree Cut",
          dendroLabels = FALSE,
          hang = 0.03,
          addGuide = TRUE,
          guideHang = 0.05,
          main = "Bicluster dendrogram and module colors"
        )
        dev.off()
        
        # Also show on console
                WGCNA::plotDendroAndColors(
                  bcTree,
                  dynamicColors,
                  "Dynamic Tree Cut",
                  dendroLabels = FALSE,
                  hang = 0.03,
                  addGuide = TRUE,
                  guideHang = 0.05,
                  main = "Bicluster dendrogram and module colors"
                )
                
      })
    } else {
      warning("Cannot plot: Length of dynamicColors does not match number of dendrogram leaves.")
    }
  }

  return(list(
    all_modules = module_bcs,
    multiomic_modules = module_bcs_multiomics,
    module_assignment = bc_module_assignment
  ))
}
