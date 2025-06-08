#' Graphical Lasso Algorithm
#'
#' This function estimates a sparse inverse covariance matrix (precision matrix) 
#' using the graphical lasso algorithm with an L1-penalized covariance estimation approach
#' by calling the fast C++ implementation from the glassoFast package.
#'
#' @param X A numeric matrix where rows represent samples and columns represent features.
#' @param rho Regularization parameter controlling sparsity.
#' @param print.detail Logical flag to print processing details (default: TRUE).
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
#' @importFrom stats cov
#' @importFrom utils capture.output
#' @importFrom glassoFast glassoFast
#' @export
graphicalLasso <- function(X, rho, print.detail = TRUE) {
  # compute the (regularized) sample covariance
  if (nrow(X) > 1) {
    Xs <- scale(X, center = TRUE, scale = TRUE)
    S  <- stats::cov(Xs)
  } else {
    # tiny ridge so glassoFast is happy
    S <- diag(ncol(X)) * 1e-4
  }
  
  # call the fast C++ implementation
  #    this will error if glassoFast is missing, but
  #    since we put it in Imports, it should always be there.
  out <- glassoFast::glassoFast(S, rho, thr = 1e-5)
  
  # return the estimated precision and covariance matrices
  list(Theta = out$wi, W = out$w)
}