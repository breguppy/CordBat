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


# ------------------
# soft threshold                                                                      Not used anymore
# ------------------
soft <- function(x, lambda){
  s <- sign(x) * max(abs(x) - lambda, 0)
  return(s)
}


# -------------------
# coordinate descent for graphical lasso sub-problem                                  Not used anymore
# -------------------
CDfgL <- function(V, beta_i, u, rho, maxIter = 200, print.detail){
  p_1 <- ncol(V)
  
  # initialize
  beta.new <- rep(0, p_1)
  finished <- rep(FALSE, p_1)
  eps <- 1e-4
  times0 <- 0
  times1 <- 0
  
  iter <- 0
  
  while(TRUE){
    iter <- iter + 1
    beta.old <- beta_i
    for (j in seq_len(p_1)) {
      df <- V %*% beta_i - u
      x <- beta_i[j] - df[j] / V[j, j]
      beta_i[j] <- soft(x, rho / V[j, j])
    }
    beta.new <- beta_i
    
    zeroIdx <- (beta_i == 0)
    if (any(zeroIdx)) {
      if (all(beta.new[zeroIdx] == beta.old[zeroIdx])) {
        times0 <- times0 + 1
      } else {
        times0 <- 0
      }
      if (times0 >= 5) {
        finished[zeroIdx] <- TRUE
      }
    }
    
    if (any(!zeroIdx)) {
      rel_change <- abs((beta.new[!zeroIdx] - beta.old[!zeroIdx]) / beta.old[!zeroIdx])
      if (all(rel_change < eps)) {
        times1 <- times1 + 1
      } else {
        times1 <- 0
      }
      if (times1 >= 3) {
        finished[!zeroIdx] <- TRUE
      }
    }
    
    if (all(finished)) break
    if (iter >= maxIter) {
      if (print.detail) message("CDfgL reached max iteration; breaking.")
      break
    }
  }
  return(beta.new)
}

# -------------------------------------------------------------
# find best parameters ksi and gamma using EBIC
# -------------------------------------------------------------
old_findBestPara <- function(X0.glist, X1.glist, penal.rho, eps, print.detail) {
  G <- length(X0.glist)
  p <- ncol(X0.glist[[1]])
  
  N0_gvec <- rep(0, G)
  N1_gvec <- rep(0, G)
  N_gvec <- rep(0, G)
  for (g in 1:G) {
    N0_gvec[g] <- nrow(X0.glist[[g]])
    N1_gvec[g] <- nrow(X1.glist[[g]])
    N_gvec[g] <- N0_gvec[g] + N1_gvec[g]
  }
  
  Sel.ksi <- 0
  Sel.gamma <- 0
  
  ksi_candidates <- c(1, 0.5, 0.3, 0.1)
  gamma_candidates <- c(1, 0.5, 0.3, 0.1)
  
  MinAvedist <- Inf
  
  for (ksi in ksi_candidates) {
    for (gamma in gamma_candidates) {
      Allpara <- BEgLasso(X0.glist, X1.glist, penal.rho, ksi, gamma, eps, print.detail)
      X1.cor.glist <- Allpara$X1.cor
      Theta.list <- Allpara$Theta
      
      r <- 0.5
      
      ebic.gvec <- rep(0, G)
      for (i in 1:G) {
        X.gi <- rbind(X0.glist[[i]], X1.cor.glist[[i]])
        X.gi.sca <- scale(X.gi, center = TRUE, scale = TRUE)
        S_i <- cov(X.gi.sca)
        Theta_i <- Theta.list[[i]]
        E.num.gi <- sum(Theta_i[upper.tri(Theta_i, diag = FALSE)] != 0)
        
        # Check determinant of Theta_i before taking its log
        detTheta <- det(Theta_i)
        if (is.na(detTheta) || detTheta <= 0) {
          ebic.gvec[i] <- Inf
        } else {
          ebic.gvec[i] <- - N_gvec[i] * (log(detTheta) - tr(S_i %*% Theta_i)) +
            E.num.gi * log(N_gvec[i]) + 4 * E.num.gi * r * log(p)
        }
      }
      
      ebic <- sum(ebic.gvec)
      
      if (!is.na(ebic) && ebic < MinAvedist) {
        MinAvedist <- ebic
        Sel.ksi <- ksi
        Sel.gamma <- gamma
      }
    }
  }
  
  penterm <- list(penal.ksi = Sel.ksi,
                  penal.gamma = Sel.gamma,
                  MinAvedist = MinAvedist)
  return(penterm)
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
