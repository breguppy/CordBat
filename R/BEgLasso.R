#' Batch Effect Correction using Group Graphical Lasso + glassoFast
#'
#' This function applies a batch effect correction method using the 
#' Group Graphical Lasso (GGM) approach, but replaces the inner coordinate-
#' descent with a fast Fortran-accelerated glasso implementation (glassoFast).
#'
#' @param X0.glist A list of reference-batch matrices (each is samples × features).
#' @param X1.glist A list of matrices to correct (each is samples × features).
#' @param penal.rho Regularization parameter for glasso.
#' @param penal.ksi Regularization parameter for coefficient update.
#' @param penal.gamma Additional penalty parameter for coefficient update.
#' @param eps Convergence threshold (default 1e-4).
#' @param print.detail Logical flag to print iteration details.
#'
#' @return A list with elements:
#'   \item{Theta}{List of precision matrices for each group.}
#'   \item{X1.cor}{List of corrected matrices for group 1.}
#'   \item{coef.a}{Scaling coefficients.}
#'   \item{coef.b}{Offset coefficients.}
#'
#' @importFrom stats cov
#' @importFrom glassoFast glassoFast
BEgLasso <- function(X0.glist,
                     X1.glist,
                     penal.rho,
                     penal.ksi,
                     penal.gamma,
                     eps = 1e-4,
                     print.detail = TRUE) {
  if (!requireNamespace("glassoFast", quietly = TRUE)) {
    stop("Please install the 'glassoFast' package to use BEgLasso().")
  }
  
  G <- length(X0.glist)
  p <- ncol(X0.glist[[1]])
  N1_gvec <- sapply(X1.glist, nrow)
  
  # Compute initial scaling/offset coefficients
  X0.m.gmat <- do.call(cbind, lapply(X0.glist, colMeans))
  X1.m.gmat <- do.call(cbind, lapply(X1.glist, colMeans))
  coef.a <- rowMeans(X0.m.gmat / X1.m.gmat)
  coef.b <- numeric(p)
  
  iter <- 0
  repeat {
    iter <- iter + 1
    
    # 1) Correct group-1 data using current coefs
    X1.cor.glist <- lapply(seq_len(G), function(g) {
      A <- diag(coef.a)
      Bmat <- matrix(coef.b, nrow = N1_gvec[g], ncol = p, byrow = TRUE)
      X1.glist[[g]] %*% A + Bmat
    })
    
    # 2) Compute covariances of combined data
    S.list <- lapply(seq_len(G), function(g) {
      Z <- rbind(X0.glist[[g]], X1.cor.glist[[g]])
      cov(scale(Z, center = TRUE, scale = TRUE))
    })
    
    # 3) Fast graphical lasso via glassoFast
    Theta.list <- lapply(S.list, function(Sg) {
      fit <- glassoFast::glassoFast(Sg, penal.rho, trace = print.detail)
      fit$wi  # precision matrix
    })
    
    # 4) Update coefficients
    upd <- update.CorrectCoef(
      X0.glist, X1.glist,
      Theta.list, coef.a, coef.b,
      penal.ksi, penal.gamma,
      print.detail
    )
    new.a <- upd$coef.a
    new.b <- upd$coef.b
    
    # 5) Check convergence
    da <- max(abs(new.a - coef.a))
    db <- max(abs(new.b - coef.b))
    coef.a <- new.a
    coef.b <- new.b
    
    if (print.detail) {
      message(sprintf("Iteration %d: |Δa|=%.5g, |Δb|=%.5g",
                      iter, da, db))
    }
    if (max(da, db) < eps) break
  }
  
  # Ensure non-negative scaling
  coef.a[coef.a < 0] <- eps
  
  return(list(
    Theta  = Theta.list,
    X1.cor = X1.cor.glist,
    coef.a = coef.a,
    coef.b = coef.b
  ))
}
