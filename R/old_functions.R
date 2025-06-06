#------------------
# Old implementation of Graphical Lasso
#------------------
#' @importFrom stats cov
#' @importFrom utils capture.output
old_graphicalLasso <- function(X, rho, print.detail) {
  N <- nrow(X)
  p <- ncol(X)
  
  if (N > 1) { # When there are multiple samples, center and scale normally. 
    # Standardize data: center and scale
    X <- scale(X, center = TRUE, scale = TRUE) 
    # Compute covariance matrix
    S <- cov(X) 
  } else { # With one sample, standardization isn’t possible since the sample variance is undefined.
    # Define S as a zero matrix (since there is no variability)  
    S <- matrix(0, p, p) 
  }
  
  # Initialize variables
  Theta <- matrix(0, p, p)
  W <- S + rho * diag(p)  # Regularized covariance matrix
  B <- matrix(0, p - 1, p)
  
  threshold <- 1e-5  # Convergence threshold
  
  # Iterative update
  repeat {
    W_old <- W  # Store previous iteration
    
    if (print.detail) {
      # Capture all output generated during the loop execution
      detailOutput <- capture.output(
        for (i in seq_len(p)) {
          idx <- setdiff(seq_len(p), i)  # Exclude current index
          
          W_11 <- W[idx, idx]
          s_12 <- S[idx, i]
          
          # Update B using penalized regression
          B[, i] <- CDfgL(W_11, B[, i], s_12, rho, print.detail = print.detail)
          W[idx, i] <- W_11 %*% B[, i]
          W[i, idx] <- W[idx, i]
        }, type = "message")
      if(length(detailOutput) > 0) {
        # Print the entire captured output as one message
        message(paste(unique(detailOutput)))
      }
    } else {
      for (i in seq_len(p)) {
        idx <- setdiff(seq_len(p), i)  # Exclude current index
        
        W_11 <- W[idx, idx]
        s_12 <- S[idx, i]
        
        # Update B using penalized regression
        B[, i] <- CDfgL(W_11, B[, i], s_12, rho, print.detail = print.detail)
        W[idx, i] <- W_11 %*% B[, i]
        W[i, idx] <- W[idx, i]
      }
    }
    
    # Convergence check
    dW <- W - W_old
    S_ndiag <- S[upper.tri(S, diag = FALSE)]
    
    # Avoid division by zero in convergence metric
    if (mean(abs(dW)) / mean(abs(S[upper.tri(S, diag = FALSE)] + 1e-6)) < threshold) break
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

#------------------
# Old implemetation of BEgLasso (with old_update_CorrectCoeff)
#------------------
#' @importFrom utils capture.output
old_BEgLasso <- function(X0.glist, 
                         X1.glist, 
                         penal.rho, 
                         penal.ksi,
                         penal.gamma, 
                         eps, 
                         print.detail) {
  
  G <- length(X0.glist)
  p <- ncol(X0.glist[[1]])
  
  # Precompute sample sizes
  N0_gvec <- sapply(X0.glist, nrow)
  N1_gvec <- sapply(X1.glist, nrow)
  
  # Compute mean matrices
  X0.m.gmat <- do.call(cbind, lapply(X0.glist, colMeans))
  X1.m.gmat <- do.call(cbind, lapply(X1.glist, colMeans))
  
  # Initialize variables
  Theta.list <- vector("list", G)
  B.list <- vector("list", G)
  coef.as <- X0.m.gmat / X1.m.gmat
  coef.a <- rowMeans(coef.as)
  coef.b <- numeric(p)
  
  coef.A <- diag(coef.a)
  
  # Initialize Theta and B for each group
  for (g in seq_len(G)) {
    Theta.list[[g]] <- matrix(0, p, p)
    B.list[[g]] <- matrix(0, p - 1, p)
  }
  
  # Correct X1 using initial coefficients
  X1.cor.glist <- lapply(seq_len(G), function(g) {
    coef.B.gi <- matrix(rep(coef.b, each = N1_gvec[g]), nrow = N1_gvec[g])
    X1.glist[[g]] %*% coef.A + coef.B.gi
  })
  
  # Compute initial covariance matrices
  W.list <- lapply(seq_len(G), function(g) {
    X.gi <- rbind(X0.glist[[g]], X1.cor.glist[[g]])
    S0_gi <- cov(scale(X.gi, center = TRUE, scale = TRUE))
    S0_gi + penal.rho * diag(1, p)
  })
  
  # Iteration variables
  finished.gmat <- matrix(FALSE, (p - 1) * p / 2, G)
  times0.gvec <- integer(G)
  times1.gvec <- integer(G)
  
  # Iterative correction
  repeat {
    W_old.list <- W.list
    
    S.list <- lapply(seq_len(G), function(g) {
      coef.B.gi <- matrix(rep(coef.b, each = N1_gvec[g]), nrow = N1_gvec[g])
      X1.gi.cor <- X1.glist[[g]] %*% coef.A + coef.B.gi
      X.gi <- rbind(X0.glist[[g]], X1.gi.cor)
      cov(scale(X.gi, center = TRUE, scale = TRUE))
    })
    
    # Update B and W using graphical lasso approach
    if (print.detail) {
      # Capture all output generated during the loop execution
      detailOutput <- capture.output(
        for (j in seq_len(p)) {
          idx <- setdiff(seq_len(p), j)
          
          for (g in seq_len(G)) {
            Wi_11 <- W.list[[g]][idx, idx]
            si_12 <- S.list[[g]][idx, j]
            
            B.list[[g]][, j] <- CDfgL(Wi_11, B.list[[g]][, j], si_12, penal.rho, 
                                      print.detail = print.detail)
            W.list[[g]][idx, j] <- Wi_11 %*% B.list[[g]][, j]
            W.list[[g]][j, idx] <- W.list[[g]][idx, j]
          }
        }
        , type = "message")
      if(length(detailOutput) > 0) {
        # Print the entire captured output as one message
        message(paste(unique(detailOutput)))
      }
    } else {
      for (j in seq_len(p)) {
        idx <- setdiff(seq_len(p), j)
        
        for (g in seq_len(G)) {
          Wi_11 <- W.list[[g]][idx, idx]
          si_12 <- S.list[[g]][idx, j]
          
          B.list[[g]][, j] <- CDfgL(Wi_11, B.list[[g]][, j], si_12, penal.rho, 
                                    print.detail = print.detail)
          W.list[[g]][idx, j] <- Wi_11 %*% B.list[[g]][, j]
          W.list[[g]][j, idx] <- W.list[[g]][idx, j]
        }
      }
    }
    
    
    # Update Theta
    for (g in seq_len(G)) {
      for (j in seq_len(p)) {
        idx <- setdiff(seq_len(p), j)
        Theta.list[[g]][j, j] <- W.list[[g]][j, j] - t(W.list[[g]][idx, j]) %*% B.list[[g]][, j]
        Theta.list[[g]][idx, j] <- -B.list[[g]][, j] * Theta.list[[g]][j, j]
        Theta.list[[g]][j, idx] <- Theta.list[[g]][idx, j]
      }
    }
    
    # Check convergence
    W_new.ndiag <- lapply(seq_len(G), function(g) W.list[[g]][upper.tri(W.list[[g]], diag = FALSE)])
    W_old.ndiag <- lapply(seq_len(G), function(g) W_old.list[[g]][upper.tri(W_old.list[[g]], diag = FALSE)])
    
    zeroIdx.list <- lapply(W_old.ndiag, function(w) w < 1e-5)
    
    for (g in seq_len(G)) {
      zeroIdx_g <- zeroIdx.list[[g]]
      if (any(zeroIdx_g)) {
        if (all(W_new.ndiag[[g]][zeroIdx_g] < 1e-5)) {
          times0.gvec[g] <- times0.gvec[g] + 1
        } else {
          times0.gvec[g] <- 0
        }
        if (times0.gvec[g] >= 5) {
          finished.gmat[zeroIdx_g, g] <- TRUE
        }
      }
      
      if (any(!zeroIdx_g)) {
        dW_g <- max(abs((W_new.ndiag[[g]][!zeroIdx_g] - W_old.ndiag[[g]][!zeroIdx_g]) 
                        / W_old.ndiag[[g]][!zeroIdx_g]))
        
        if (dW_g < eps) {
          times1.gvec[g] <- times1.gvec[g] + 1
        } else {
          times1.gvec[g] <- 0
        }
        
        if (times1.gvec[g] >= 3) {
          finished.gmat[!zeroIdx_g, g] <- TRUE
        }
      }
    }
    
    if (all(as.vector(finished.gmat))) break
    
    # Update coefficients
    coef.update <- update_CorrectCoef(X0.glist, X1.glist, Theta.list, coef.a, 
                                      coef.b, penal.ksi, penal.gamma, print.detail)
    coef.a <- coef.update$coef.a
    coef.b <- coef.update$coef.b
  }
  
  # Ensure non-negative coef.a
  coef.a[coef.a < 0] <- 0.05
  
  ### This part is new ###
  # Re-compute the final X1.cor.glist
  coef.A <- diag(coef.a)
  X1.cor.glist <- lapply(seq_len(G), function(g) {
    coef.B.gi <- matrix(rep(coef.b, each = N1_gvec[g]), nrow = N1_gvec[g])
    X1.glist[[g]] %*% coef.A + coef.B.gi
  })
  
  list(Theta = Theta.list, X1.cor = X1.cor.glist, coef.a = coef.a, coef.b = coef.b)
}


# ------------------
# soft threshold
# ------------------
soft <- function(x, lambda){
  s <- sign(x) * pmax(abs(x) - lambda, 0)
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
# old update correction coefficients a and b
# -------------------------------------------------------------
update_CorrectCoef <- function(X0.glist, X1.glist, Theta.list, 
                               a.i, b.i, penal.ksi, penal.gamma, print.detail) {
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
  
  a.o <- a.i
  b.o <- b.i
  
  for (j in seq_len(p)) {
    
    # update a[j]
    tmp1.gvec <- rep(0, G)
    tmp2.gvec <- rep(0, G)
    for (g in seq_len(G)) {
      A <- diag(a.o)
      B_gi <- matrix(rep(b.o, each = N1_gvec[g]), nrow = N1_gvec[g])
      
      X1.gi.cor <- X1.glist[[g]] %*% A + B_gi
      X.gi <- rbind(X0.glist[[g]], X1.gi.cor)
      
      X.gi.sca <- scale(X.gi, center = TRUE, scale = TRUE)
      X.gi.sca.attr <- attributes(X.gi.sca)
      Mu_g <- X.gi.sca.attr$`scaled:center`
      Sigma_g <- X.gi.sca.attr$`scaled:scale`
      
      # Ensure Sigma_g[j] is non-zero and not NA
      if (is.na(Sigma_g[j]) || abs(Sigma_g[j]) < 1e-6) {
        Sigma_g[j] <- 1e-6
        if (print.detail) warning(sprintf("Sigma_g[%d] was zero or NA for group
                                          %d in a-update; replaced with 1e-6", j, g))
      }
      
      # Use only the rows corresponding to group 1 in the scaled data
      if ((N0_gvec[g] + 1) > N_gvec[g]) {
        stop("Not enough rows in the combined data for group 1 extraction")
      }
      Y_g <- X.gi.sca[(N0_gvec[g] + 1):N_gvec[g], ]
      
      # Replace the j-th column: adjust b.o[j] by the scaling factors
      if(is.null(nrow(Y_g))) {
        Y_g[j] <- (rep(b.o[j], N1_gvec[g]) - Mu_g[j]) / Sigma_g[j]
      }else {
        Y_g[, j] <- (rep(b.o[j], N1_gvec[g]) - Mu_g[j]) / Sigma_g[j]
      }
      Z_g <- Y_g
      
      tmp <- as.vector(Z_g %*% Theta.list[[g]][, j])
      tmp1.gvec[g] <- 2 / (N_gvec[g] * Sigma_g[j]) * sum(X1.glist[[g]][, j] * tmp, na.rm = TRUE)
      tmp2.gvec[g] <- 2 * Theta.list[[g]][j, j] / (N_gvec[g] * Sigma_g[j]^2) * sum(X1.glist[[g]][, j]^2, na.rm = TRUE)
    }
    denominator_a <- sum(tmp2.gvec, na.rm = TRUE)
    if (abs(denominator_a) < 1e-6) {
      if (print.detail) warning(sprintf("Denominator for updating a[%d] is zero;
                                        keeping previous value", j))
      a.o[j] <- a.o[j]
    } else {
      a.o[j] <- 1 + soft(- sum(tmp1.gvec, na.rm = TRUE) - sum(tmp2.gvec, na.rm = TRUE), penal.ksi) / denominator_a
    }
    
    # update b[j]
    tmp3.gvec <- rep(0, G)
    tmp4.gvec <- rep(0, G)
    for (g in seq_len(G)) {
      A <- diag(a.o)
      B_gi <- matrix(rep(b.o, each = N1_gvec[g]), nrow = N1_gvec[g])
      X1.gi.cor <- X1.glist[[g]] %*% A + B_gi
      X.gi <- rbind(X0.glist[[g]], X1.gi.cor)
      
      X.gi.sca <- scale(X.gi, center = TRUE, scale = TRUE)
      X.gi.sca.attr <- attributes(X.gi.sca)
      Mu_g <- X.gi.sca.attr$`scaled:center`
      Sigma_g <- X.gi.sca.attr$`scaled:scale`
      
      if (is.na(Sigma_g[j]) || abs(Sigma_g[j]) < 1e-6) {
        Sigma_g[j] <- 1e-6
        if (print.detail) warning(sprintf("Sigma_g[%d] was zero or NA for group 
                                          %d in b-update; replaced with 1e-6", j, g))
      }
      Y_g <- X.gi.sca[(N0_gvec[g] + 1):N_gvec[g], ]
      
      # Update j-th column with the current a.o[j] and adjust by scaling
      if(is.null(nrow(Y_g))) {
        Y_g[j] <- (X1.glist[[g]][, j] * a.o[j] - Mu_g[j]) / Sigma_g[j]
      } else {
        Y_g[, j] <- (X1.glist[[g]][, j] * a.o[j] - Mu_g[j]) / Sigma_g[j]
      }
      Z_g <- Y_g
      
      tmp <- as.vector(Z_g %*% Theta.list[[g]][, j])
      tmp3.gvec[g] <- 2 / (N_gvec[g] * Sigma_g[j]) * sum(tmp, na.rm = TRUE)
      tmp4.gvec[g] <- 2 * N1_gvec[g] * Theta.list[[g]][j, j] / (N_gvec[g] * Sigma_g[j]^2)
    }
    denominator_b <- sum(tmp4.gvec, na.rm = TRUE)
    if (abs(denominator_b) < 1e-6) {
      if (print.detail) warning(sprintf("Denominator for updating b[%d] is 
                                        zero; keeping previous value", j))
      b.o[j] <- b.o[j]
    } else {
      b.o[j] <- soft(- sum(tmp3.gvec, na.rm = TRUE), penal.gamma) / denominator_b
    }
  }
  
  coef.out <- list(coef.a = a.o, coef.b = b.o)
  return(coef.out)
}


# -------------------------------------------------------------
# find best parameters for penalty selection via CV+BIC (grid search)
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
old_selrho.useCVBIC <- function(X, print.detail = TRUE) {
  N <- nrow(X)
  fold <- selfoldforCV(N)
  CVset.size <- N / fold
  
  rhos <- seq(from = 0.1, to = 0.9, by = 0.1)
  CVerr1 <- matrix(0, fold, length(rhos))
  CVerr2 <- rep(0, length(rhos))
  
  for (r in seq_along(rhos)) {
    rho <- rhos[r]
    if (fold != 1) {
      for (i in 1:fold) {
        start.index <- (i-1) * CVset.size + 1
        end.index <- i * CVset.size
        X.cv <- X[start.index:end.index, ]
        X.tr <- X[-(start.index:end.index), ]
        
        c.mat <- graphicalLasso(X.tr, rho, print.detail = print.detail)
        Theta <- c.mat$Theta
        
        X.cv.sca <- scale(X.cv, center = TRUE, scale = TRUE)
        S.cv <- cov(X.cv.sca)
        
        k <- sum(Theta[upper.tri(Theta, diag = FALSE)] != 0)
        CVerr1[i, r] <- k * log(CVset.size) - CVset.size * (log(det(Theta)) - tr(S.cv %*% Theta))
      }
    } else {
      c.mat <- graphicalLasso(X, rho, print.detail = print.detail)
      Theta <- c.mat$Theta
      # Added to handle matrix of 1 row
      if (nrow(X) > 1) { 
        X.sca <- scale(X, center = TRUE, scale = TRUE) 
        S <- cov(X.sca) 
      } else { 
        S <- matrix(0, ncol(X), ncol(X)) 
      }
      
      k <- sum(Theta[upper.tri(Theta, diag = FALSE)] != 0)
      CVerr2[r] <- k * log(N) - 2 * (log(det(Theta)) - tr(S %*% Theta))
    }
  }
  
  if (fold != 1) {
    CVerr1 <- colMeans(CVerr1)
    MinCVerr <- min(CVerr1)
    rho.cv <- rhos[CVerr1 == MinCVerr]
  } else {
    MinCVerr <- min(CVerr2)
    rho.cv <- rhos[CVerr2 == MinCVerr]
  }
  
  if (print.detail) {
    message('CVBIC: select rho =', rho.cv, '\n')
  }
  
  return(c(rho.cv, MinCVerr))
}
