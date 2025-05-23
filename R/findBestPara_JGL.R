#' Joint‐Graphical‐Lasso CordBat wrapper with explicit EBIC surface
#'
#' @param X0.glist     List of reference-batch matrices (each n×p).
#' @param X1.glist     List of to-be-corrected matrices (each m×p).
#' @param lambda1.seq  Numeric vector of glasso penalties (λ1).
#' @param lambda2.seq  Numeric vector of fusion penalties (λ2).
#' @param r             EBIC gamma parameter (default 0.5).
#'
#' @return A list with:
#'   \item{Theta}{List of G selected precision matrices.}
#'   \item{best}{Named numeric of length 2 (lambda1, lambda2) with custom subsetting.}
#'   \item{EBIC}{Matrix of summed EBIC values (length(lambda1.seq) × length(lambda2.seq)).}
#'
#' @importFrom JGL JGL
#' @export
findBestPara_JGL <- function(X0.glist,
                             X1.glist,
                             lambda1.seq,
                             lambda2.seq,
                             r = 0.5) {
  if (!requireNamespace("JGL", quietly = TRUE)) {
    stop("Please install the 'JGL' package to use findBestPara_JGL().")
  }
  G <- length(X0.glist)
  p <- ncol(X0.glist[[1]])
  
  # Prepare combined data per group
  Y.list <- lapply(seq_len(G), function(g) {
    rbind(X0.glist[[g]], X1.glist[[g]])
  })
  
  # Initialize EBIC surface
  L1 <- length(lambda1.seq)
  L2 <- length(lambda2.seq)
  ebic.surface <- matrix(NA_real_, nrow = L1, ncol = L2,
                         dimnames = list(
                           paste0("lambda1=", lambda1.seq),
                           paste0("lambda2=", lambda2.seq)
                         ))
  
  best.ebic   <- Inf
  best.p1     <- NA_real_
  best.p2     <- NA_real_
  best.Thetas <- NULL
  
  # Grid search
  for (i1 in seq_along(lambda1.seq)) {
    lam1 <- lambda1.seq[i1]
    for (i2 in seq_along(lambda2.seq)) {
      lam2 <- lambda2.seq[i2]
      
      out <- JGL::JGL(
        Y                  = Y.list,
        penalty            = "group",
        lambda1            = lam1,
        lambda2            = lam2,
        weights = "equal",
        return.whole.theta = TRUE
      )
      Thetas <- out$theta  # list-of-G precision matrices
      
      # Compute EBIC for this pair
      tot <- 0
      for (g in seq_len(G)) {
        Zg <- Y.list[[g]]
        ng <- nrow(Zg)
        Sg <- cov(scale(Zg, center = TRUE, scale = TRUE))
        Θ  <- Thetas[[g]]
        E  <- sum(Θ[upper.tri(Θ)] != 0)
        tot <- tot + (
          - ng * (log(det(Θ)) - sum(Sg * Θ))
          + E * log(ng)
          + 4 * E * r * log(p)
        )
      }
      ebic.surface[i1, i2] <- tot
      
      if (tot < best.ebic) {
        best.ebic   <- tot
        best.p1     <- lam1
        best.p2     <- lam2
        best.Thetas <- Thetas
      }
    }
  }
  
  # Build named best vector with custom class
  best <- c(lambda1 = best.p1, lambda2 = best.p2)
  class(best) <- "cordbat_best"
  
  list(
    Theta = best.Thetas,
    best  = best,
    EBIC  = ebic.surface
  )
}

#' Subsetting method to drop names when extracting single penalty value
#'
#' Prevent name mismatch when doing best["lambda1"] in tests
#'
#' @export
`[.cordbat_best` <- function(x, i, ...) {
  val <- NextMethod()
  if (is.numeric(val) && length(val) == 1) {
    names(val) <- NULL
  }
  val
}
