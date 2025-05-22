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
# soft threshold
# ------------------
soft <- function(x, lambda){
  s <- sign(x) * max(abs(x) - lambda, 0)
  return(s)
}


# -------------------
# coordinate descent for graphical lasso sub-problem
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
# find best parameters for penalty selection via CV+BIC (grid search)
# -------------------------------------------------------------
findBestPara <- function(X0.glist, X1.glist, penal.rho, eps, print.detail) {
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
# EBIC select rho via huge::huge + huge::huge.select
# -------------------------------------------------------------
#' @importFrom stats cov
#' @importFrom lava tr
selrho.useCVBIC <- function(X,
                            print.detail = TRUE,
                            nlambda      = 50,
                            gamma        = 0) {
  if (!requireNamespace("huge", quietly = TRUE)) {
    stop("Please install the 'huge' package to use EBIC selection.")
  }
  
  # 1) Fit the full glasso path
  out <- huge::huge(
    X,
    method  = "glasso",
    nlambda = nlambda,
    verbose = print.detail
  )
  
  # 2) EBIC selection (note EBIC.gamma and positional first arg)
  sel <- huge::huge.select(
    out,
    criterion   = "ebic",
    ebic.gamma  = gamma,
    verbose     = print.detail
  )
  
  # 3) Extract the chosen lambda
  rho.ebic <- sel$opt.lambda
  
  # 4) Grab the EBIC‐vector from the correct slot (‐> sel$ebic)
  ebic.vec <- if (!is.null(sel$ebic)) {
    sel$ebic
  } else {
    stop("Unexpected huge.select output: no sel$ebic component found")
  }
  ebic.val <- ebic.vec[ sel$opt.index ]
  
  # 5) Message & return
  if (print.detail) {
    message(sprintf(
      "EBIC: select rho = %s  (EBIC = %.2f)",
      rho.ebic, ebic.val
    ))
  }
  
  return(c(rho.ebic, ebic.val))
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
