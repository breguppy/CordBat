#' StARS via huge::huge + huge::huge.select
#'
#' @param X             numeric matrix (samples × features)
#' @param b             subsample size for stability (integer)
#' @param M             number of subsampling replications (integer)
#' @param seed          optional RNG seed for reproducibility
#' @param beta          instability threshold (default 0.05)
#' @param print.detail  logical; if TRUE, prints progress
#' @return numeric vector c(Sel.rho, D_var)
#' @importFrom huge huge huge.select
StARS_huge <- function(X, b, M,
                       seed = NULL,
                       beta    = 0.05,
                       print.detail = TRUE) {
  if (!requireNamespace("huge", quietly = TRUE)) {
    stop("Please install the 'huge' package to use StARS_huge().")
  }
  if (!is.null(seed)) set.seed(seed)
  n <- nrow(X)
  # sub-sample ratio must be in (0,1]
  subsamp_ratio <- b / n
  
  candidate_grids <- list(
    seq(0.9, 0.5, by = -0.1),
    seq(0.5, 0.1, by = -0.1),
    seq(0.1, 0.05, by = -0.01),
    seq(0.05, 0.01, by = -0.01)
  )
  for (grid in candidate_grids) {
    # Fit glasso
    out <- huge::huge(X,
                      method = "glasso",
                      lambda = grid,
                      cov.output = FALSE,
                      verbose = print.detail)
    # StARS selection
    sel <- huge::huge.select(out,
                             criterion = "stars",
                             stars.thresh = beta,
                             stars.subsample.ratio = subsamp_ratio,
                             rep.num = M,
                             verbose = print.detail)
  
    # pull out the chosen rho (i.e. lambda) and its instability measure
    Sel.rho <- sel$opt.lambda
    D_var   <- sel$variability[sel$opt.index]
  
    if (D_var <= beta) {
      if (print.detail) {
        message("StARS selected rho = ", Sel.rho,
                " with instability = ", round(D_var, 4))
      }
      return(c(Sel.rho = Sel.rho, D_var = D_var))
    }
  }
  
  # If no rho satisfied the beta threshold, return the last attempted
  if (print.detail) {
    warning("No rho met the instability threshold; returning best found.")
  }
  return(c(Sel.rho = Sel.rho, D_var = D_var))
}
