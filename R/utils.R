# -------------------------------------------------------------
# select a proper fold number for CV
# -------------------------------------------------------------
selfoldforCV <- function(N){
  foldstosel <- 2:9
  Num.quo <- N %/% foldstosel
  Num.rem <- N %% foldstosel
  
  # if no fold meets requirements, then change N
  selIdx <- (Num.rem == 0 & Num.quo >= 10)
  folds <- foldstosel[selIdx]
  folds.num <- length(folds)
  
  if (folds.num != 0) {
    d5 <- abs(rep(5, folds.num) - folds)
    fold <- folds[d5 == min(d5)]
    if (length(fold) > 1) {
      if (N > 300) {
        fold <- max(fold)
      } else {
        fold <- min(fold)
      }
    }
  } else {
    fold <- 1
  }
  
  return(fold)
}


#----------------------------------------------------
# findBestPara: cached & vectorized EBIC grid search
#----------------------------------------------------
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
  
  # pick minimal EBIC with tie-breaking: largest ksi, then largest gamma
  best_val  <- min(ebic_vals)
  best_idxs <- which(ebic_vals == best_val)
  best_df   <- params[best_idxs, , drop=FALSE]
  # order by increasing ksi then gamma
  tie_order <- order(-best_df$ksi, -best_df$gamma)
  chosen    <- best_idxs[tie_order[1]]
  
  list(
    penal.ksi   = params$ksi[chosen],
    penal.gamma = params$gamma[chosen],
    MinAvedist  = ebic_vals[chosen]
  )
}


# -------------------------------------------------------------
# CV + BIC select rho
# -------------------------------------------------------------
#' @importFrom stats cov
#' @importFrom lava tr
selrho.useCVBIC <- function(X, print.detail = TRUE) {
  N <- nrow(X)
  if (N <= 1) {
    if (print.detail) message("selrho.useCVBIC(): need at least 2 rows to do CV+BIC selection.")
    return(c(rho.cv = 0.1, MinCVerr = NA_real_))
  }
  fold <- selfoldforCV(N)
  CVset.size <- N / fold
  
  rhos <- seq(from = 0.1, to = 0.9, by = 0.1)
  R <- length(rhos)
  
  CVerr1 <- matrix(0, nrow = fold, ncol = R)
  CVerr2 <- numeric(R)
  
  # mask for the “upper triangle” of a p×p matrix
  p <- ncol(X)
  mask <- upper.tri(matrix(NA, p, p))
  
  if (fold > 1) {
    for (r in seq_len(R)) {
      rho <- rhos[r]
      
      for (i in seq_len(fold)) {
        # split out CV fold i
        start.index <- (i - 1) * CVset.size + 1
        end.index   <- i * CVset.size
        
        X.cv <- X[start.index:end.index, , drop = FALSE]
        X.tr <- X[-(start.index:end.index), , drop = FALSE]
        
        # fit glasso on the training fold
        Theta <- graphicalLasso(X.tr, rho, print.detail = print.detail)$Theta
        
        # scale+cov of the held‐out fold
        X.cv.sca <- scale(X.cv, center = TRUE, scale = TRUE)
        S.cv     <- cov(X.cv.sca)
        
        # count edges in the upper triangle of Theta
        k <- sum(Theta[mask] != 0)
        
        # exactly the same BIC‐style error for that fold:
        CVerr1[i, r] <- k * log(CVset.size) -
          CVset.size * (log(det(Theta)) - tr(S.cv %*% Theta))
      }
    }
    
    # average over folds, pick the rho that gives the smallest mean
    CVerr.mean <- colMeans(CVerr1)
    MinCVerr  <- min(CVerr.mean)
    idx_all <- which(CVerr.mean == MinCVerr)
    i_max   <- max(idx_all)          # index of the largest‐rho among the ties
    rho.cv  <- rhos[i_max]
    
  } else {
    # fold = 1
    # scale the entire X once
    X.sca <- scale(X, center = TRUE, scale = TRUE)
    S     <- cov(X.sca)
    
    for (r in seq_len(R)) {
      Theta <- graphicalLasso(X, rhos[r])$Theta
      k <- sum(Theta[mask] != 0)
      CVerr2[r] <- k * log(N) -
        2 * (log(det(Theta)) - tr(S %*% Theta))
    }
    
    MinCVerr <- min(CVerr2)
    idx_all <- which(CVerr2 == MinCVerr)
    i_max   <- max(idx_all)          # index of the largest‐rho among the ties
    rho.cv  <- rhos[i_max]
  }
  
  if (print.detail) {
    message(
      "CV + BIC selects rho = ", rho.cv,
      " with MinCVerr = ", round(MinCVerr, 4)
    )
  }
  
  return(c(rho.cv, MinCVerr))
}


# -------------------------------------------------------------
# Delete outliers in data 
# -------------------------------------------------------------
#' @importFrom stats sd
#' @importFrom mixOmics pca
DelOutlier <- function(X) {
  pca.dat <- pca(X, ncomp = 3, center = TRUE, scale = TRUE)
  pca.dat.varX <- pca.dat$variates$X
  delsampIdx <- c()
  for (i in 1:3) {
    pc.i <- pca.dat.varX[, i]
    pc.i.m <- mean(pc.i)
    pc.i.sd <- sd(pc.i)
    pc.i.min <- pc.i.m - 3 * pc.i.sd
    pc.i.max <- pc.i.m + 3 * pc.i.sd
    delsampIdx <- c(delsampIdx, which(pc.i < pc.i.min | pc.i > pc.i.max))
  }
  delsampIdx <- unique(delsampIdx)
  
  if (length(delsampIdx) == 0) {
    X.out <- X
  } else {
    X.out <- X[-delsampIdx, ]
  }
  
  Del.result <- list(delsampIdx = delsampIdx, 
                     X.out = X.out)
  return(Del.result)
}


# -------------------------------------------------------------
# Impute outliers in data
# -------------------------------------------------------------
#' @importFrom stats sd
#' @importFrom DMwR2 knnImputation
ImputeOutlier <- function(X) {
  p <- ncol(X)
  X.out <- X
  for (i in 1:p) {
    dat.i <- X[, i]  # Fixed use column i instead of p
    dat.i.m <- mean(dat.i)
    dat.i.sd <- sd(dat.i)
    dat.i.max <- dat.i.m + 4 * dat.i.sd
    dat.i.min <- dat.i.m - 4 * dat.i.sd
    dat.i[dat.i < dat.i.min | dat.i > dat.i.max] <- NA
    X.out[, i] <- dat.i  # assign to column i
  }
  
  na.num <- sum(is.na(X.out))
  if (na.num != 0) {
    X.out <- as.data.frame(X.out)
    X.out <- knnImputation(X.out)
  }
  
  X.out <- as.matrix(X.out)
  
  return(X.out)
}
