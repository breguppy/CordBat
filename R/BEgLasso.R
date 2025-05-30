#' Batch Effect Correction using Group Graphical Lasso (updated)
#'
#' This function applies a batch effect correction method using the
#' Group Graphical Lasso (GGM) approach, leveraging the fast C++ backend
#' via the `graphicalLasso()` helper for covariance and precision updates.
#'
#' @param X0.glist A list of matrices where each matrix corresponds to a batch
#'                  in group 0 (reference batch).
#' @param X1.glist A list of matrices where each matrix corresponds to a batch
#'                  in group 1 (batch to be corrected).
#' @param penal.rho Regularization parameter for the graphical lasso.
#' @param penal.ksi Regularization parameter for coefficient update.
#' @param penal.gamma Additional penalty parameter for coefficient update.
#' @param eps Convergence threshold for stopping criteria.
#' @param print.detail Logical flag to print processing details (default: TRUE).
#'
#' @return A list containing:
#'   \itemize{
#'     \item \code{Theta} - List of precision matrices for each group.
#'     \item \code{X1.cor} - List of corrected matrices for group 1.
#'     \item \code{coef.a} - Corrected scaling coefficients.
#'     \item \code{coef.b} - Corrected offset coefficients.
#'   }
#'
#' @importFrom stats cov
#' @importFrom utils capture.output
#' @export
BEgLasso <- function(
    X0.glist,
    X1.glist,
    penal.rho,
    penal.ksi,
    penal.gamma,
    eps,
    print.detail = TRUE
) {
  G <- length(X0.glist)
  p <- ncol(X0.glist[[1]])
  
  # Precompute sample sizes
  N1_gvec <- sapply(X1.glist, nrow)
  
  # Initial coefficient estimates from group means
  X0.means <- do.call(cbind, lapply(X0.glist, colMeans))
  X1.means <- do.call(cbind, lapply(X1.glist, colMeans))
  coef.a <- rowMeans(X0.means / X1.means)
  coef.b <- numeric(p)
  
  # Helper to compute corrected X1 for each group
  make_X1_cor <- function(a, b) {
    lapply(seq_len(G), function(g) {
      Bgi <- matrix(rep(b, each = N1_gvec[g]), nrow = N1_gvec[g])
      X1.glist[[g]] %*% diag(a) + Bgi
    })
  }
  
  # Initial corrected data
  X1.cor.glist <- make_X1_cor(coef.a, coef.b)
  
  # Initial graphical-lasso fit for W & Theta
  W.list <- vector("list", G)
  Theta.list <- vector("list", G)
  for (g in seq_len(G)) {
    Xg <- rbind(X0.glist[[g]], X1.cor.glist[[g]])
    gl <- graphicalLasso(Xg, rho = penal.rho, print.detail = print.detail)
    W.list[[g]]     <- gl$W
    Theta.list[[g]] <- gl$Theta
  }
  
  # Iteration bookkeeping
  finished.gmat <- matrix(FALSE, (p - 1) * p / 2, G)
  times0.gvec <- integer(G)
  times1.gvec <- integer(G)
  
  # Iterative correction loop
  repeat {
    W_old.list <- W.list
    
    # Compute scaled covariances with current correction
    S.list <- lapply(seq_len(G), function(g) {
      Xg <- rbind(
        X0.glist[[g]],
        X1.glist[[g]] %*% diag(coef.a) +
          matrix(rep(coef.b, each = N1_gvec[g]), nrow = N1_gvec[g])
      )
      cov(scale(Xg, center = TRUE, scale = TRUE))
    })
    
    # Graphical-lasso update per group
    for (g in seq_len(G)) {
      # Build combined data for group g
      Xg <- rbind(
        X0.glist[[g]],
        X1.glist[[g]] %*% diag(coef.a) +
          matrix(rep(coef.b, each = N1_gvec[g]), nrow = N1_gvec[g])
      )
      gl <- graphicalLasso(Xg, rho = penal.rho, print.detail = print.detail)
      W.list[[g]]     <- gl$W
      Theta.list[[g]] <- gl$Theta
    }
    
    # Check convergence on off-diagonals
    W_new_nd <- lapply(W.list, function(W) W[upper.tri(W, diag = FALSE)])
    W_old_nd <- lapply(W_old.list, function(W) W[upper.tri(W, diag = FALSE)])
    zeroIdx <- lapply(W_old_nd, function(v) v < eps)
    
    for (g in seq_len(G)) {
      zi <- zeroIdx[[g]]
      # track entries that have effectively shrunk to zero
      if (any(zi)) {
        if (all(W_new_nd[[g]][zi] < eps)) {
          times0.gvec[g] <- times0.gvec[g] + 1
        } else {
          times0.gvec[g] <- 0
        }
        if (times0.gvec[g] >= 5) finished.gmat[zi, g] <- TRUE
      }
      nz <- !zi
      if (any(nz)) {
        dW <- max(abs((W_new_nd[[g]][nz] - W_old_nd[[g]][nz]) /
                        W_old_nd[[g]][nz]))
        if (dW < eps) {
          times1.gvec[g] <- times1.gvec[g] + 1
        } else {
          times1.gvec[g] <- 0
        }
        if (times1.gvec[g] >= 3) finished.gmat[nz, g] <- TRUE
      }
    }
    
    if (all(as.vector(finished.gmat))) break
    
    # Update a & b via fast C++ routine
    upd <- update.CorrectCoef(
      X0.glist, X1.glist, Theta.list,
      coef.a, coef.b,
      penal.ksi, penal.gamma,
      print.detail
    )
    coef.a <- upd$coef.a
    coef.b <- upd$coef.b
  }
  
  # Clamp negative a's
  coef.a[coef.a < 0] <- 0.05
  
  # Final corrected X1
  X1.cor.glist <- make_X1_cor(coef.a, coef.b)
  
  list(
    Theta = Theta.list,
    X1.cor = X1.cor.glist,
    coef.a = coef.a,
    coef.b = coef.b
  )
}
