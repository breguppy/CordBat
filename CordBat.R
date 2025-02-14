#' Batch Effect Correction using CordBat
#'
#' This function applies batch effect correction (BEC) using a community detection-based 
#' approach (CordBat). It handles outlier removal, imputation, and correction across batches.
#'
#' @param X A numeric matrix of log-transformed expression data (samples x features).
#' @param batch A vector of batch identifiers corresponding to rows in `X`.
#' @param group An optional vector of group identifiers corresponding to rows in `X` (default: NULL).
#' @param grouping Logical indicating whether to use provided groupings (default: FALSE).
#' @param ref.batch The reference batch identifier.
#' @param eps Convergence threshold for batch effect correction (default: 1e-5).
#' @param print.detail Logical flag to print processing details (default: TRUE).
#' 
#' @return A list containing:
#' \itemize{
#'   \item \code{batch.level} - Unique batch levels.
#'   \item \code{delsampIdx} - Indices of deleted samples (outliers).
#'   \item \code{batch.new} - Updated batch assignments after outlier removal.
#'   \item \code{group.new} - Updated group assignments after outlier removal.
#'   \item \code{X.delout} - Data matrix after outlier removal.
#'   \item \code{X.cor} - Final corrected data matrix.
#'   \item \code{X.cor.1} - Corrected data excluding outliers.
#'   \item \code{X.cor.withQC} - Corrected data including QC samples (if present).
#'   \item \code{Xcor.para} - List of batch-specific correction parameters.
#' }
#' 
#' @examples
#' set.seed(123)
#' X <- matrix(rnorm(1000), 100, 10)
#' batch <- rep(1:5, each = 20)
#' res <- CordBat(X, batch, ref.batch = 1, print.detail = FALSE)
#' str(res)
#'
CordBat <- function(X, batch, group = NULL, grouping = FALSE, ref.batch, 
                    eps = 1e-5, print.detail = TRUE) {
  
  # If group is not provided, assign all to a single group
  if (is.null(group)) {
    containQC <- NA
    group <- rep(1, nrow(X))
  } else {
    containQC <- any(group == "QC")
    
    if (containQC) {
      X.init <- X
      batch.init <- batch
      group.init <- group
      
      # Separate QC samples
      X_QC <- X[group == "QC", ]
      QC_batch <- batch[group == "QC"]
      
      # Remove QC samples for processing
      X <- X[group != "QC", ]
      batch <- batch[group != "QC"]
      group <- group[group != "QC"]
    }
    
    if (!grouping) group <- rep(1, nrow(X))
  }
  
  X <- as.matrix(X)
  p <- ncol(X)
  
  # Convert batch and group to factors
  batch.f <- factor(batch)
  batch.levels <- levels(batch.f)
  batch.num <- length(batch.levels)
  
  group.f <- factor(group)
  group.levels <- levels(group.f)
  group.num <- length(group.levels)
  
  # Outlier Removal and Imputation
  delsampIdx <- integer()
  X.delout <- X
  
  for (i in seq_len(batch.num)) {
    bati.idx <- which(batch == batch.levels[i])
    X.bati <- DelOutlier(X[bati.idx, ])
    delsamp.bati <- X.bati$delsampIdx
    dat.bati <- ImputeOutlier(X.bati$X.out)
    
    if (length(delsamp.bati) > 0) {
      bati.delinitIdx <- bati.idx[delsamp.bati]
      delsampIdx <- c(delsampIdx, bati.delinitIdx)
      X.delout[bati.idx[-delsamp.bati], ] <- dat.bati
    } else {
      X.delout[bati.idx, ] <- dat.bati
    }
  }
  
  if (length(delsampIdx) > 0) {
    X.nodel <- X.delout
    X.delout <- X.delout[-delsampIdx, ]
    batch <- batch[-delsampIdx]
    group <- group[-delsampIdx]
  } else {
    X.nodel <- X
  }
  
  batch.f <- factor(batch)
  group.f <- factor(group)
  
  # Initialize correction parameters
  Theta.list <- vector("list", group.num)
  for (i in seq_len(group.num)) Theta.list[[i]] <- matrix(0, p, p)
  
  a <- rep(1, p)
  b <- rep(0, p)
  para <- list(Theta = Theta.list, coef.a = a, coef.b = b)
  
  Xcor.para <- vector("list", batch.num)
  for (i in seq_len(batch.num)) Xcor.para[[i]] <- para
  
  # Initialize corrected matrices
  X.cor <- matrix(0, nrow(X), ncol(X))
  X.cor[batch == ref.batch, ] <- X[batch == ref.batch, ]
  
  X.cor.1 <- matrix(0, nrow(X.delout), ncol(X.delout))
  X.cor.1[batch == ref.batch, ] <- X.delout[batch == ref.batch, ]
  
  X.cor.withQC <- if (containQC) {
    X.tmp <- matrix(0, nrow(X.init), ncol(X.init))
    X.tmp[batch.init == ref.batch, ] <- X.init[batch.init == ref.batch, ]
    X.tmp
  } else NULL
  
  # Community detection
  Xb0.mat <- X.delout[batch == ref.batch, ]
  COM <- getAllCom(Xb0.mat)
  
  if (print.detail) cat("Community detection:", length(COM), "communities\n")
  
  # Batch effect correction
  for (i in seq_along(COM)) {
    metID <- COM[[i]]
    Xb0.COMi.glist <- lapply(seq_len(group.num), function(g) {
      X.delout[which(batch == ref.batch & group == group.levels[g]), metID]
    })
    
    rhos <- sapply(seq_len(group.num), function(g) {
      if (length(metID) > 5) {
        StARS(Xb0.COMi.glist[[g]], round(0.7 * nrow(Xb0.COMi.glist[[g]])), 100, print.detail)[1]
      } else {
        selrho.useCVBIC(Xb0.COMi.glist[[g]], print.detail)[1]
      }
    })
    
    rho <- mean(rhos)
    if (print.detail) cat('Set rho =', rho, '\n')
    
    for (k in seq_len(batch.num)) {
      if (batch.levels[k] == ref.batch) next
      
      Xb1.Batk.COMi.glist <- lapply(seq_len(group.num), function(g) {
        X.delout[which(batch == batch.levels[k] & group == group.levels[g]), metID]
      })
      
      penterm <- findBestPara(Xb0.COMi.glist, Xb1.Batk.COMi.glist, rho, eps)
      
      para.out <- BEgLasso(Xb0.COMi.glist, Xb1.Batk.COMi.glist, rho, penterm$penal.ksi, penterm$penal.gamma, eps)
      
      for (g in seq_len(group.num)) {
        X.cor.1[which(batch == batch.levels[k] & group == group.levels[g]), metID] <- para.out$X1.cor[[g]]
      }
    }
    
    if (print.detail) cat('Finished correction of community', i, '\n')
  }
  
  list(batch.level = batch.levels, delsampIdx = delsampIdx, batch.new = batch, 
       group.new = group, X.delout = X.delout, X.cor = X.cor, X.cor.1 = X.cor.1, 
       X.cor.withQC = X.cor.withQC, Xcor.para = Xcor.para)
}
