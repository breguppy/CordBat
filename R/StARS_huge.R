#' StARS via huge::huge + huge::huge.select
#'
#' @param X        Numeric data matrix (n samples × p features)
#' @param b        Subsample size for each stability run
#' @param M        Number of subsampling iterations (rep.num)
#' @param nlambda  How many λ’s to fit on the path (default 50)
#' @param beta     Stability threshold (stars.thresh), default 0.05
#' @param verbose  Print progress if TRUE
#' @return Named vector c(Sel.rho, D_var)
StARS_huge <- function(X, b, M,
                       nlambda = 50,
                       beta    = 0.05,
                       verbose = TRUE) {
  if (!requireNamespace("huge", quietly = TRUE)) {
    stop("Please install the 'huge' package to use StARS_huge().")
  }
  
  n <- nrow(X)
  # translate your fixed subsample size to a ratio in (0,1]
  subsamp_ratio <- b / n
  
  # 1) fit a full glasso path in C/C++
  out <- huge::huge(X,
                    method   = "glasso",
                    nlambda  = nlambda,
                    verbose  = verbose)
  
  # 2) perform StARS selection
  sel <- huge::huge.select(out,
                           criterion             = "stars",
                           stars.thresh          = beta,
                           stars.subsample.ratio = subsamp_ratio,
                           rep.num               = M,
                           verbose               = verbose)
  
  # 3) pull out the chosen λ (i.e. ρ) and its instability measure
  Sel.rho <- sel$opt.lambda
  D_var   <- sel$variability[sel$opt.index]
  
  if (verbose) {
    message("StARS selected rho = ", Sel.rho,
            " with instability = ", round(D_var, 4))
  }
  
  return(c(Sel.rho = Sel.rho, D_var = D_var))
}
