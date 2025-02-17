#' Graphical Lasso Algorithm
#'
#' This function estimates a sparse inverse covariance matrix (precision matrix) 
#' using the graphical lasso algorithm with an L1-penalized covariance estimation approach.
#'
#' @param X A numeric matrix where rows represent samples and columns represent features.
#' @param rho Regularization parameter controlling sparsity.
#' 
#' @return A list containing:
#' \itemize{
#'   \item \code{Theta} - The estimated precision (inverse covariance) matrix.
#'   \item \code{W} - The estimated covariance matrix.
#' }
#' 
#' @examples
#' set.seed(123)
#' X <- matrix(rnorm(100), 10, 10)
#' result <- graphicalLasso(X, rho = 0.1)
#' str(result)
#' 
graphicalLasso <- function(X, rho) {
  N <- nrow(X)
  p <- ncol(X)
  
  # Standardize data: center and scale
  X <- scale(X, center = TRUE, scale = TRUE)
  
  # Compute covariance matrix
  S <- cov(X)
  
  # Initialize variables
  Theta <- matrix(0, p, p)
  W <- S + rho * diag(p)  # Regularized covariance matrix
  B <- matrix(0, p - 1, p)
  
  threshold <- 1e-5  # Convergence threshold
  
  # Iterative update
  repeat {
    W_old <- W  # Store previous iteration
    
    for (i in seq_len(p)) {
      idx <- setdiff(seq_len(p), i)  # Exclude current index
      
      W_11 <- W[idx, idx]
      s_12 <- S[idx, i]
      
      # Update B using penalized regression
      B[, i] <- CDfgL(W_11, B[, i], s_12, rho)
      W[idx, i] <- W_11 %*% B[, i]
      W[i, idx] <- W[idx, i]
    }
    
    # Convergence check
    dW <- W - W_old
    S_ndiag <- S[upper.tri(S, diag = FALSE)]
    
    if (mean(abs(dW)) / mean(abs(S_ndiag)) < threshold) break
  }
  
  # Compute precision matrix Theta
  for (i in seq_len(p)) {
    idx <- setdiff(seq_len(p), i)
    Theta[i, i] <- 1 / (W[i, i] - t(W[idx, i]) %*% B[, i])
    Theta[idx, i] <- -B[, i] * Theta[i, i]
    Theta[i, idx] <- Theta[idx, i]
  }
  
  # Return results
  list(Theta = Theta, W = W)
}
