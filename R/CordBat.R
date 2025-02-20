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
#'   \item \code{X.cor} - Fully corrected data matrix (using outlier-free imputed data).
#'   \item \code{X.cor.1} - Corrected data matrix (using outlier-removed data).
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
#' @export
CordBat <- function(X, 
                    batch, 
                    group = NULL, 
                    grouping = FALSE, 
                    ref.batch, 
                    eps = 1e-5, 
                    print.detail = TRUE) {
  
  # If group is not provided, assign all samples to a single group
  if (is.null(group)) {
    containQC <- FALSE
    group <- rep(1, nrow(X))
  } else {
    containQC <- any(group == "QC")
    if (containQC) {
      X.init <- X
      batch.init <- batch
      group.init <- group
      
      # Separate QC samples
      X_QC <- X[group == "QC", , drop = FALSE]
      QC_batch <- batch[group == "QC"]
      
      # Remove QC samples for processing
      X <- X[group != "QC", , drop = FALSE]
      batch <- batch[group != "QC"]
      group <- group[group != "QC"]
    }
    if (!grouping) group <- rep(1, nrow(X))
  }
  
  X <- as.matrix(X)
  p <- ncol(X)
  
  # Force batch and group to be character vectors for comparisons
  batch <- as.character(batch)
  group <- as.character(group)
  
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
    cur.batch <- batch.levels[i]
    bati.idx <- which(batch == cur.batch)
    if (length(bati.idx) == 1) next
    X.bati <- DelOutlier(X[bati.idx, , drop = FALSE])
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
    #X.delout <- X.delout[-delsampIdx, , drop = FALSE]
    #batch <- batch[-delsampIdx]
    #group <- group[-delsampIdx]
  } else {
    X.nodel <- X
  }
  
  # Refresh factors after outlier removal
  batch.f <- factor(batch)
  batch.levels <- levels(batch.f)
  batch.num <- length(batch.levels)
  group.f <- factor(group)
  group.levels <- levels(group.f)
  group.num <- length(group.levels)
  
  # Initialize correction parameters
  Theta.list <- vector("list", group.num)
  for (i in seq_len(group.num)) Theta.list[[i]] <- matrix(0, p, p)
  a <- rep(1, p)
  b <- rep(0, p)
  para <- list(Theta = Theta.list, coef.a = a, coef.b = b)
  
  Xcor.para <- vector("list", batch.num)
  for (i in seq_len(batch.num)) Xcor.para[[i]] <- para
  
  # Initialize corrected matrices (dimensions based on outlier-removed data)
  #n_del <- nrow(X.delout)
  n_del <- nrow(X)
  X.cor   <- matrix(0, n_del, p)
  X.cor.1 <- matrix(0, n_del, p)
  
  ref.batch_char <- as.character(ref.batch)
  ref.idx <- which(batch == ref.batch_char)
  X.cor[ref.idx, ]   <- X.delout[ref.idx, ]
  X.cor.1[ref.idx, ] <- X.delout[ref.idx, ]
  
  X.cor.withQC <- NULL
  if (containQC) {
    X.cor.withQC <- matrix(0, nrow(X.init), p)
    ref.idx.init <- which(as.character(batch.init) == as.character(ref.batch))
    X.cor.withQC[ref.idx.init, ] <- X.init[ref.idx.init, ]
  }
  
  # === Community detection on reference batch data ===
  ref_batch_idx <- which(batch == ref.batch_char)
  Xb0.mat <- X.delout[ref_batch_idx, , drop = FALSE]
  COM <- getAllCom(Xb0.mat)
  if (print.detail) message("Community detection: ", length(COM), " communities")
  
  # === Batch effect correction for each community ===
  for (i in seq_along(COM)) {
    metID <- COM[[i]]
    
    # Build list for reference batch communities by group (skip groups with no data)
    Xb0.COMi.glist <- vector("list", group.num)
    valid_groups_ref <- logical(group.num)
    for (g in seq_len(group.num)) {
      grp_label <- group.levels[g]
      idx <- which(batch == ref.batch_char & group == grp_label)
      if (length(idx) > 0) {
        Xb0.COMi.glist[[g]] <- X.delout[idx, metID, drop = FALSE]
        valid_groups_ref[g] <- TRUE
      } else {
        valid_groups_ref[g] <- FALSE
      }
    }
    # Keep only valid groups
    Xb0.COMi.glist <- Xb0.COMi.glist[valid_groups_ref]
    if (length(Xb0.COMi.glist) == 0) next  # skip community if no valid group exists
    
    # Determine rho from the reference groups
    rhos <- sapply(Xb0.COMi.glist, function(mat) {
      if (nrow(mat) > 5) {
        StARS(mat, round(0.7 * nrow(mat)), 100, print.detail)[1]
      } else {
        selrho.useCVBIC(mat, print.detail)[1]
      }
    })
    rho <- mean(rhos, na.rm = TRUE)
    if (print.detail) message("Set rho = ", rho)
    
    # Process each non-reference batch
    for (k in seq_len(batch.num)) {
      batch_label <- batch.levels[k]
      if (batch_label == ref.batch_char) next
      if (length(which(batch == batch_label)) < 2) next
      
      # Build list for non-reference batch communities by group
      Xb1.Batk.COMi.glist <- vector("list", group.num)
      valid_groups_nonref <- logical(group.num)
      for (g in seq_len(group.num)) {
        grp_label <- group.levels[g]
        idx <- which(batch == batch_label & group == grp_label)
        if (length(idx) > 0) {
          Xb1.Batk.COMi.glist[[g]] <- X.delout[idx, metID, drop = FALSE]
          valid_groups_nonref[g] <- TRUE
        } else {
          valid_groups_nonref[g] <- FALSE
        }
      }
      Xb1.Batk.COMi.glist <- Xb1.Batk.COMi.glist[valid_groups_nonref]
      if (length(Xb1.Batk.COMi.glist) == 0) next  # skip if no valid groups
      
      # Use the valid reference and non-reference groups in downstream functions
      penterm <- findBestPara(Xb0.COMi.glist, Xb1.Batk.COMi.glist, rho, eps)
      para.out <- BEgLasso(Xb0.COMi.glist, Xb1.Batk.COMi.glist, rho, 
                           penterm$penal.ksi, penterm$penal.gamma, eps)
      
      # Update corrected data (X.cor.1) for this non-reference batch
      # We update for each valid group; note that the order in Xb1.Batk.COMi.glist now corresponds 
      # to the subset of groups with data (which we can retrieve from valid_groups_nonref)
      valid_grp_indices <- which(valid_groups_nonref)
      for (j in seq_along(valid_grp_indices)) {
        g_idx <- valid_grp_indices[j]
        idx <- which(batch == batch_label & group == group.levels[g_idx])
        X.cor.1[idx, metID] <- para.out$X1.cor[[j]]
      }
      
      # Update fully corrected data (X.cor) using outlier-free data (X.nodel)
      idx_nodel <- which(batch == batch_label)
      if (length(idx_nodel) > 0) {
        Xb1.nodel <- X.nodel[idx_nodel, , drop = FALSE]
        N1 <- nrow(Xb1.nodel)
        coef.A <- diag(para.out$coef.a)
        coef.B <- matrix(rep(para.out$coef.b, each = N1), nrow = N1)
        X.cor[idx_nodel, metID] <- Xb1.nodel[, metID] %*% coef.A + coef.B
      }
    }
    
    if (print.detail) message("Finished correction of community ", i)
  }
  
  return(list(batch.level = batch.levels, 
              delsampIdx = delsampIdx, 
              batch.new = batch, 
              group.new = group, 
              X.delout = X.delout, 
              X.cor = X.cor, 
              X.cor.1 = X.cor.1, 
              X.cor.withQC = X.cor.withQC, 
              Xcor.para = Xcor.para))
}
