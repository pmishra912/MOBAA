#' find_ensemble_bc
#'
#' @description Wrapper function for mosbi biclustering with user-selected methods
#'
#' @param data_matrix Omic data matrix with features as rows and samples as columns
#' @param methods A character vector of biclustering methods to use. Options: "fabia", "isa", "plaid", "qubic"
#'
#' @return List of ensemble biclusters.
#' @examples
#' find_ensemble_bc(data_matrix, methods = c("fabia", "qubic"))
#' @export
find_ensemble_bc <- function(data_matrix, methods = c("fabia", "isa", "plaid", "qubic"),quiet = TRUE,seed = NULL) {
    
    if (!is.null(seed)) {
        RNGkind(kind = "L'Ecuyer-CMRG")
        set.seed(seed)
      }
    
    if (!requireNamespace("mosbi", quietly = TRUE)) {
        stop("The 'mosbi' package is required but not installed.")
      }

  all_bics <- list()

  run_safe <- function(expr) {
      if (quiet) {
        suppressMessages(suppressWarnings(capture.output(res <- try(expr, silent = TRUE))))
      } else {
        res <- try(expr, silent = FALSE)
      }
      return(res)
    }
  
  if ("fabia" %in% methods) {
    bcfb <- run_safe(mosbi::run_fabia(data_matrix))
    all_bics <- c(all_bics, bcfb)
  }

  if ("isa" %in% methods) {
    bcisa <- run_safe(mosbi::run_isa(data_matrix))
    all_bics <- c(all_bics, bcisa)
  }

  if ("plaid" %in% methods) {
    bcplaid <- run_safe(mosbi::run_plaid(data_matrix))
    all_bics <- c(all_bics, bcplaid)
  }

  if ("qubic" %in% methods) {
    bcqubic <- run_safe(mosbi::run_qubic(data_matrix))
    all_bics <- c(all_bics, bcqubic)
  }

  # Check if any biclusters were generated
  if (length(all_bics) == 0) {
    stop("No valid biclustering methods specified or all methods failed.")
  }

  # Suppress potential plot outputs if quiet = TRUE
    if (quiet) {
      tmp_plot <- tempfile(fileext = ".png")
      png(tmp_plot)  # divert plots
    }
    
  # Continue with network and ensemble
  bic_net <- mosbi::bicluster_network(
    all_bics,
    data_matrix,
    n_randomizations = 5,
    MARGIN = "both",
    metric = 4,
    n_steps = 1000
  )

  coms <- mosbi::get_louvain_communities(
    bic_net,
    min_size = 3,
    bics = all_bics
  )

  ensemble_bicluster_list <- run_safe(mosbi::ensemble_biclusters(
    coms,
    all_bics,
    data_matrix,
    row_threshold = 0.1,
    col_threshold = 0.1
  ))

  if (quiet) {
      dev.off()  # close diverted plot
      unlink(tmp_plot)  # delete temp plot file
    }
  
  return(ensemble_bicluster_list)
}
