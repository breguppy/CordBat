# ----------------------------------------
# select a proper fold number for CV
# --------------------------------------
selfoldforCV <- function(N){
  foldstosel <- c(2:9)
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
# coordinate descent 
# -------------------
CDfgL <- function(V, beta_i, u, rho){
  p_1 <- ncol(V)
  
  # initialize
  beta.new <- rep(0, p_1)
  finished <- rep(FALSE, p_1)
  eps <- 1e-4
  times0 <- 0
  times1 <- 0
  
  while(TRUE){
    beta.old <- beta_i
    for (j in c(1: p_1)) {
      df <- V %*% beta_i - u
      x <- beta_i[j] - df[j] / V[j, j]
      beta_i[j] <- soft(x, rho / V[j, j])
    }
    beta.new <- beta_i
    
    zeroIdx <- (beta_i == 0)
    if (any(zeroIdx)) {
      if (all(beta.new[zeroIdx] == beta.old[zeroIdx])) {
        times0 <- times0 + 1
      } else{
        times0 <- 0
      }
      if (times0 >= 5) {
        finished[zeroIdx] <- TRUE
      }
    }
    
    if (any(!zeroIdx)) {
      if (all(abs((beta.new[!zeroIdx] - beta.old[!zeroIdx])
                  / beta.old[!zeroIdx]) < eps)) {
        times1 <- times1 + 1
      } else{
        times1 <- 0
      }
      if (times1 >= 3) {
        finished[!zeroIdx] <- TRUE
      }
    }
    if (all(finished)) {
      break
    }
  }
  return(beta.new)
}

# --------------------------------------------
# CV + BIC select rho
# --------------------------------------------
selrho.useCVBIC <- function(X, print.detail = T) {
  N <- nrow(X)
  fold <- selfoldforCV(N)
  CVset.size <- N / fold
  
  rhos <- seq(from = 0.1, to = 0.9, by = 0.1)
  CVerr1 <- matrix(0, fold, length(rhos))
  CVerr2 <- rep(0, length(rhos))
  
  for (r in c(1: length(rhos))) {
    rho <- rhos[r]
    if (fold != 1) {
      for (i in c(1: fold)) {
        start.index <- (i-1) * CVset.size + 1
        end.index <- i * CVset.size
        X.cv <- X[c(start.index: end.index), ]
        X.tr <- X
        X.tr <- X.tr[-c(start.index: end.index), ]
        
        c.mat <- graphicalLasso(X.tr, rho)
        Theta <- c.mat$Theta
        
        # compute error for CV set
        X.cv.sca <- scale(X.cv, center = TRUE, scale = TRUE)
        S.cv <- cov(X.cv.sca)
        
        k <- sum(Theta[upper.tri(Theta, diag = FALSE)] != 0)
        CVerr1[i, r] <- k * log(CVset.size) - CVset.size * (log(det(Theta)) - tr(S.cv %*% Theta))
      }
    } else {
      c.mat <- graphicalLasso(X, rho)
      Theta <- c.mat$Theta
      
      # compute error for CV set
      X.sca <- scale(X, center = TRUE, scale = TRUE)
      S <- cov(X.sca)
      
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
    cat('CVBIC: select rho =', rho.cv, '\n')
  }
  
  return(c(rho.cv, MinCVerr))
}

# -------------------------
# Delete outliers in data 
# --------------------------
DelOutlier <- function(X) {
  pca.dat <- pca(X, ncomp = 3, center = TRUE, scale = TRUE)
  pca.dat.varX <- pca.dat$variates$X
  delsampIdx <- c()
  for (i in c(1: 3)) {
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

# --------------------------
# Impute outliers in data
# -------------------------
ImputeOutlier <- function(X) {
  p <- ncol(X)
  X.out <- X
  for (i in c(1: p)) {
    dat.i <- X[, p]
    dat.i.m <- mean(dat.i)
    dat.i.sd <- sd(dat.i)
    dat.i.max <- dat.i.m + 3 * dat.i.sd
    dat.i.min <- dat.i.m - 3 * dat.i.sd
    dat.i[dat.i < dat.i.min | dat.i > dat.i.max] <- NA
    
    X.out[, p] <- dat.i
  }
  
  na.num <- sum(is.na(X.out))
  if (na.num != 0) {
    X.out <- as.data.frame(X.out)
    X.out <- knnImputation(X.out)
  }
  
  X.out <- as.matrix(X.out)
  
  return(X.out)
}

# -------------------
# get reference batch
# ---------------------
getRefbat <- function(Data, batch) {
  # cumulative RSD
  batch <- as.factor(batch)
  batch.num <- nlevels(batch)
  batch.lev <- levels(batch)
  
  rsd <- seq(0.01, 1, 0.01)
  rsd.cumfreq <- data.frame(batch = factor(rep(batch.lev, each = length(rsd)), 
                                           levels = c(1: batch.num)), 
                            rsd = rep(rsd, batch.num), 
                            cumfreq = rep(0, length(rsd) * batch.num))
  
  p <- ncol(Data)
  for (i in c(1: batch.num)) {
    dat.i <- Data[batch == batch.lev[i], ]
    m.val <- apply(dat.i, MARGIN = 2, FUN = mean)
    sd.val <- apply(dat.i, MARGIN = 2, FUN = sd)
    rsd.val <- sd.val / m.val
    cumfreq.i <- rep(0, length(rsd))
    
    for (k in c(1: length(rsd))) {
      cumfreq.i[k] <- sum(rsd.val <= rsd[k]) / p
    }
    rsd.cumfreq$cumfreq[rsd.cumfreq$batch == batch.lev[i]] <- cumfreq.i
  }
  
  # compute the area under the curve
  AUCFC_rsd <- rep(0, batch.num) # area under the cumulative frequency curve of RSD for each batch
  delta_x <- 0.01
  for (i in c(1: batch.num)) {
    cumfreq.bati <- rsd.cumfreq[rsd.cumfreq$batch == batch.lev[i], ]$cumfreq
    n <- length(cumfreq.bati)
    
    y.1 <- cumfreq.bati[1]
    y.n <- cumfreq.bati[n]
    y.is <- cumfreq.bati[-c(1, n)]
    
    AUCFC_rsd[i] <- delta_x * (0.5 * (y.1 + y.n) + sum(y.is))
  }
  
  batch.rank <- as.numeric(batch.lev[order(AUCFC_rsd, decreasing = T)])
  AUCFC_rsd.sort <- sort(AUCFC_rsd, decreasing = T)
  
  # plot
  linesize <- 1
  cbPalette <- c("#E69F00", "#56B4E9", "#009E73", "#CC79A7", 
                 "#F0E442", "#0072B2", "#D55E00", "#999999", 
                 "#33CC99", "#990000", "#CC99FF", "#FF00FF", 
                 '#5500FF', '#66FF66', '#0088A8', '#A500CC')
  
  lab.txt <- paste0('Area under the curve: \n(decreasing order)\n')
  for (i in c(1: batch.num)) {
    lab.txt <- paste0(lab.txt, 
                      'Batch ', batch.rank[i], ': ', round(AUCFC_rsd.sort[i], 4), '\n')
  }
  
  ggplot(data = rsd.cumfreq, aes(x = rsd, y = cumfreq)) +
    geom_line(aes(color = batch), linewidth = linesize) +
    geom_vline(xintercept = 0.2, linetype = 2) +
    geom_vline(xintercept = 0.3, linetype = 2) + 
    geom_text(aes(x = 0.85, y = 0.02 * (batch.num + 2)), 
              label = lab.txt, color = "black", size = 5) + 
    xlab('RSD') + ylab('Cumulative Frequency') + xlim(0, 1) + ylim(0, 1) +
    scale_color_manual(values = cbPalette) + scale_linetype_discrete() + theme_bw() +
    theme(axis.title = element_text(size = 15), axis.text = element_text(size = 12, face = 'bold'),
          legend.title = element_text(size = 15), legend.text = element_text(size = 12),
          plot.title = element_text(size = 18, hjust = 0.5)) +
    labs(color = 'Batch', linetype = NA, title = 'RSD')
}