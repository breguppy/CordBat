#' StARS via huge::huge + huge::huge.select
#'
#' @param X             numeric matrix (samples × features)
#' @param b             subsample size for stability (integer)
#' @param M             number of subsampling replications (integer)
#' @param lambda.grid   vector of rhos to try (e.g. seq(0.9,0.1,-0.1) or seq(0.1,0.01,-0.01))
#' @param seed          optional RNG seed for reproducibility
#' @param beta          instability threshold (default 0.05)
#' @param print.detail  logical; if TRUE, prints progress
#' @return numeric vector c(Sel.rho, D_var)
#' @importFrom huge huge huge.select
StARS_huge <- function(X, b, M,
                       lambda.grid,
                       seed = NULL,
                       beta    = 0.05,
                       print.detail = TRUE) {
  if (!requireNamespace("huge", quietly = TRUE)) {
    stop("Please install the 'huge' package to use StARS_huge().")
  }
  if (!is.null(seed)) set.seed(seed)
  n <- nrow(X)
  # translate your fixed subsample size to a ratio in (0,1]
  subsamp_ratio <- b / n
  
  # 1) fit a full glasso path in C/C++
  out <- huge::huge(X,
                    method   = "glasso",
                    lambda  = lambda.grid,
                    cov.output     = FALSE,
                    verbose        = print.detail)
  
  # 2) perform StARS selection
  sel <- huge::huge.select(out,
                           criterion             = "stars",
                           stars.thresh          = beta,
                           stars.subsample.ratio = subsamp_ratio,
                           rep.num               = M,
                           verbose               = print.detail)
  
  # 3) pull out the chosen rho (i.e. lambda) and its instability measure
  Sel.rho <- sel$opt.lambda
  D_var   <- sel$variability[sel$opt.index]
  
  if (print.detail) {
    message("StARS selected rho = ", Sel.rho,
            " with instability = ", round(D_var, 4))
  }
  
  return(c(Sel.rho = Sel.rho, D_var = D_var))
}
