#' findBestPara: cached & vectorized EBIC grid search
#'
#' Grid–search over (ksi, gamma) using BEgLasso calls,
#' then EBIC‐score each fit with a vectorized mapply.
#'
#' @param X0.glist  List of reference‐batch data matrices.
#' @param X1.glist  List of batch‐to‐correct data matrices.
#' @param penal.rho Graphical‐lasso penalty.
#' @param eps       Convergence tolerance passed to BEgLasso.
#' @param print.detail Logical; forward to BEgLasso (default: FALSE).
#'
#' @return A list with the best \code{penal.ksi}, \code{penal.gamma}, and \code{MinAvedist}.
#' @export
findBestPara <- function(X0.glist, X1.glist, penal.rho, eps, print.detail=FALSE) {
  # special-case: identical data => default (1,1)
  if (identical(X0.glist, X1.glist)) {
    return(list(
      penal.ksi   = 1,
      penal.gamma = 1,
      MinAvedist  = 0
    ))
  }
  
  G <- length(X0.glist)
  p <- ncol(X0.glist[[1]])
  r <- 0.5
  logp <- log(p)
  
  # sample sizes
  N0 <- vapply(X0.glist, nrow, integer(1))
  N1 <- vapply(X1.glist, nrow, integer(1))
  N  <- N0 + N1
  
  # parameter grid
  ksis   <- c(1, 0.5, 0.3, 0.1)
  gammas <- ksis
  params <- expand.grid(ksi=ksis, gamma=gammas, KEEP.OUT.ATTRS=FALSE)
  nGrid  <- nrow(params)
  
  # store EBIC for each grid point
  ebic_vals <- numeric(nGrid)
  
  # evaluate EBIC across grid
  for (idx in seq_len(nGrid)) {
    ksi   <- params$ksi[idx]
    gamma <- params$gamma[idx]
    
    fit   <- BEgLasso(X0.glist, X1.glist, penal.rho, ksi, gamma, eps, print.detail)
    X1c   <- fit$X1.cor
    Theta <- fit$Theta
    
    # vectorized EBIC per group
    ebic_g <- mapply(function(X0, X1cg, Th, Ni) {
      Xgi   <- rbind(X0, X1cg)
      Si    <- cov(scale(Xgi, center=TRUE, scale=TRUE))
      detTh <- det(Th)
      if (is.na(detTh) || detTh <= 0) return(Inf)
      E     <- sum(Th[upper.tri(Th)] != 0)
      -Ni * (log(detTh) - sum(Si * Th)) + E * log(Ni) + 4 * E * r * logp
    }, X0.glist, X1c, Theta, N, SIMPLIFY=TRUE)
    
    ebic_vals[idx] <- sum(ebic_g)
  }
  
  # pick minimal EBIC with tie-breaking: smallest ksi, then smallest gamma
  best_val  <- min(ebic_vals)
  best_idxs <- which(ebic_vals == best_val)
  best_df   <- params[best_idxs, , drop=FALSE]
  # order by increasing ksi then gamma
  tie_order <- order(best_df$ksi, best_df$gamma)
  chosen    <- best_idxs[tie_order[1]]
  
  list(
    penal.ksi   = params$ksi[chosen],
    penal.gamma = params$gamma[chosen],
    MinAvedist  = ebic_vals[chosen]
  )
}