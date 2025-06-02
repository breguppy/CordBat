#' CV + BIC selector via CVglasso & huge::huge.select
#'
#' @param X            numeric matrix (samples × features)
#' @param rhos         numeric vector of penalties (e.g. seq(0.1,0.9,by=0.1))
#' @param print.detail logical; if TRUE, prints a summary message
#' @param seed         optional integer seed for reproducibility
#' @return named numeric vector c(rho.sel, MinCVerr)
#' @importFrom CVglasso CVglasso
#' @importFrom huge huge huge.select
select_rho_cv_bic <- function(X,
                              rhos,
                              seed = NULL,
                              print.detail = TRUE) {
  if (!requireNamespace("CVglasso", quietly = TRUE)) {
    stop("Please install the 'CVglasso' package to use this function.")
  }
  if (!requireNamespace("huge", quietly = TRUE)) {
    stop("Please install the 'huge' package to use this function.")
  }
  # 1) determine folds
  N    <- nrow(X)
  fold <- selfoldforCV(N)   # your existing helper for choosing fold
  
  # 2) set seed if requested
  if (!is.null(seed)) set.seed(seed)
  
  # 3) CV case
  if (fold > 1) {
    cvfit <- CVglasso::CVglasso(
      X,
      K           = fold,
      lam           = rhos,
      diagonal = FALSE,
      tol         = 1e-5
    )
    rho.sel   <- cvfit$Tuning[2]
    MinCVerr  <- cvfit$MIN.error 
    
    # 4) BIC fallback
  } else {
    huge_out <- huge::huge(
      X,
      method            = "glasso",
      lambda            = rhos,
      cov.output        = FALSE,
      verbose           = print.detail
    )
    sel <- huge::huge.select(
      huge_out,
      criterion    = "ebic",
      ebic.gamma   = 0,
      verbose      = print.detail
    )
    rho.sel  <- sel$opt.lambda                 # chosen lambda
    MinCVerr <- sel$ebic[sel$opt.index] # BIC score at lambda
  }
  
  # 5) print if requested
  if (print.detail) {
    message(
      "CV + BIC selects rho = ", rho.sel,
      " with min avg CV error = ", round(MinCVerr, 4)
    )
  }
  
  # 6) return both values
  return(c(rho.sel = rho.sel, MinCVerr = MinCVerr))
}