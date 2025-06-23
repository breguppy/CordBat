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
#' @param skip.impute Logical flag to detect, delete and impute outliers (default: FALSE)
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
#' @importFrom utils capture.output
#' @export
CordBat <- function(X, 
                    batch, 
                    group = NULL, 
                    grouping = FALSE, 
                    ref.batch, 
                    eps = 1e-5, 
                    print.detail = TRUE, 
                    skip.impute = FALSE) {
  
  # If group is not provided, assign all samples to a single group
  if (is.null(group)) {
    containQC <- FALSE
    group <- rep(1, nrow(X))
  } else {
    # If group is provided, remove QC samples for processing
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
  
  batch.old <- batch
  group.old <- group
  delsampIdx <- integer()
  if(!skip.impute) {
    # Outlier Removal and Imputation
    batch.new <- batch
    group.new <- group
    X.delout <- X
    for (i in seq_len(batch.num)) {
      cur.batch <- batch.levels[i]
      bati.idx <- which(batch == cur.batch)
      if (length(bati.idx) == 1) next  # Skip batches with only 1 sample
      X.bati <- DelOutlier(X[bati.idx, , drop = FALSE]) # Determine sample outliers with PCA
      delsamp.bati <- X.bati$delsampIdx
      
      # impute outlier measurements in non-sample outliers.
      if (length(delsamp.bati) > 0) {
        dat.bati <- ImputeOutlier(X.bati$X.out)  # Perform imputation only if deletion happened
        bati.delinitIdx <- bati.idx[delsamp.bati]
        delsampIdx <- c(delsampIdx, bati.delinitIdx)
        X.delout[bati.idx[-delsamp.bati], ] <- dat.bati
      } else {
        X.delout[bati.idx, ] <- X.bati$X.out  # No deletion, use original data
      }
    }
    if (length(delsampIdx) > 0) {
      X.nodel <- X.delout                    # imputed data with sample outliers remaining untouched
      X.delout <- X.delout[-delsampIdx, ]    # imputed data with sample outliers removed
      batch.new <- batch.new[-delsampIdx]
      group.new <- group.new[-delsampIdx]
    } else {
      X.nodel <- X
    }
    # Refresh factors after outlier removal
    batch <- batch.new
    batch.f <- factor(batch)
    batch.levels <- levels(batch.f)
    batch.num <- length(batch.levels)
    group <- group.new
    group.f <- factor(group)
    group.levels <- levels(group.f)
    group.num <- length(group.levels)
  } else {
    X.nodel <- X # Skip outlier detection/imputation
    X.delout <- X
  }
  
  # Initialize correction parameters
  Theta.list <- vector("list", group.num)
  for (i in seq_len(group.num)) Theta.list[[i]] <- matrix(0, p, p)
  a <- rep(1, p)
  b <- rep(0, p)
  para <- list(Theta = Theta.list, coef.a = a, coef.b = b)
  
  Xcor.para <- vector("list", batch.num)
  for (i in seq_len(batch.num)) Xcor.para[[i]] <- para
  
  # Initialize corrected matrices 
  X.cor <- matrix(0, nrow(X.nodel), p)
  X.cor.1 <- matrix(0, nrow(X.delout), p)
  
  ref.batch_char <- as.character(ref.batch)
  ref.idx.old <- which(batch.old == ref.batch_char)
  ref.idx <- which(batch == ref.batch_char)
  X.cor[ref.idx.old, ] <- X.nodel[ref.idx.old, ]
  X.cor.1[ref.idx, ] <- X.delout[ref.idx, ]
  
  
  X.cor.withQC <- NULL
  if ((!is.na(containQC)) & containQC) {
    X.cor.withQC <- matrix(0, nrow(X.init), p)
    ref.idx.init <- which(as.character(batch.init) == as.character(ref.batch))
    X.cor.withQC[ref.idx.init, ] <- X.init[ref.idx.init, ]
  }
  
  # === Community detection on reference batch data ===
  ref_batch_idx <- which(batch == ref.batch_char)
  Xb0.mat <- X.delout[ref_batch_idx, , drop = FALSE]
  COM <- getAllCom(Xb0.mat)
  if (print.detail) message("Community detection: ", length(COM), " communities", "\n",
                            "Size: ", paste(vapply(COM, length, integer(1)), collaspe = " "), "\n")
  
  numCores <- parallel::detectCores()
  if (print.detail) {
    message("Detected ", numCores, " logical cores on this machine.")
  }
  
  cl <- parallel::makeCluster(parallel::detectCores())
  on.exit(parallel::stopCluster(cl), add = TRUE)
  
  ## double‐check the cluster size
  numWorkers <- length(cl)
  if (print.detail) {
    message("Spun up a PSOCK cluster with ", numWorkers, " workers.")
  }
  
  ## export all needed data & functions
  parallel::clusterExport(cl,
                          varlist = c("COM", "batch", "batch.old", "batch.levels", "batch.num",
                                      "group", "group.levels", "group.num",
                                      "X.delout", "X.nodel", "ref.batch_char",
                                      "eps", "print.detail"),
                          envir = environment()
  )
  parallel::clusterEvalQ(cl, { 
    library(CordBat)        # for findBestPara, BEgLasso, StARS, selrho.useCVBIC, etc.
    library(utils) 
  })
  
  # === Batch effect correction for each community ===
  comm_res <- parallel::parLapply(cl, seq_along(COM), function(i) {
    metID <- COM[[i]]
    
    # Build list for reference batch communities by group (skip groups with no data)
    X0_glist <- vector("list", group.num)
    valid0   <- logical(group.num)
    for (g in seq_len(group.num)) {
      idx <- which(batch == ref.batch_char & group == group.levels[g])
      if (length(idx) > 0) {
        X0_glist[[g]] <- X.delout[idx, metID, drop = FALSE]
        valid0[g] <- TRUE
      } else {
        valid0[g] <- FALSE
      }
    }
    # Keep only valid groups
    X0_glist <- X0_glist[valid0]
    if (length(X0_glist) == 0) return(NULL)  # skip community if no valid group exists
    
    # Determine rho from the reference groups
    pick_rho <- function(mat) {
      ## if someone passed you a 0×p matrix, bail out
      if (nrow(mat) == 0) return( NA_real_ )
      
      ## choose which routine to call
      safe_val <- tryCatch({
        if (nrow(mat) > 5) {
          StARS(mat, round(0.7 * nrow(mat)), 100, print.detail)[1]
        } else {
          selrho.useCVBIC(mat, print.detail)[1]
        }
      }, error = function(e) {
        warning("pick_rho(): caught error, returning NA: ", e$message)
        NA_real_
      })
      
      ## if they somehow gave you length-0 or NULL, coerce
      if (length(safe_val) != 1 || is.null(safe_val)) {
        safe_val <- NA_real_
      }
      
      ## emit the StARS/CVBIC messages if requested
      if (print.detail) {
        ## re-capture and print only once
        detailOutput <- capture.output({
          if (nrow(mat) > 5) {
            StARS(mat, round(0.7 * nrow(mat)), 100, print.detail)
          } else {
            selrho.useCVBIC(mat, print.detail)
          }
        }, type = "message")
        message(paste(unique(detailOutput), collapse = "\n"))
      }
      
      return(safe_val)
    }
    
    rhos <- vapply(X0_glist, pick_rho, numeric(1))
    rho  <- mean(rhos, na.rm = TRUE)
    if (print.detail) message(" Community ", i, ": rho=", rho)
    
    # Process each non-reference batch
    local_para <- vector("list", batch.num)
    for (k in seq_len(batch.num)) {
      lbl <- batch.levels[k]
      if (lbl == ref.batch_char) next
      if (length(which(batch == lbl)) < 2) next
      
      # Build list for non-reference batch communities by group
      X1_glist <- vector("list", group.num)
      valid1 <- logical(group.num)
      for (g in seq_len(group.num)) {
        grp_label <- group.levels[g]
        idx <- which(batch == lbl & group == grp_label)
        if (length(idx) > 0) {
          X1_glist[[g]] <- X.delout[idx, metID, drop = FALSE]
          valid1[g] <- TRUE
        } else {
          valid1[g] <- FALSE
        }
      }
      X1_glist <- X1_glist[valid1]
      if (length(X1_glist) == 0) next  # skip if no valid groups
      
      # Use the valid reference and non-reference groups in downstream functions
      penterm <- findBestPara(X0_glist, X1_glist, rho, eps, print.detail)
      
      if (print.detail) {
        cat('Batch ', k, ' correction begin......')
      }
      
      para.out <- BEgLasso(X0_glist, X1_glist, rho, 
                           penterm$penal.ksi, penterm$penal.gamma, eps, 
                           print.detail)
      if (print.detail) {
        cat('finshed', '\n')
        if (anyNA(para.out$coef.a) || anyNA(para.out$coef.b)) {
          stop("Got NA in correction coefficients; aborting.")
        }
      }
      
      local_para[[k]] <- list(
        metID       = metID,
        valid1      = valid1,
        coef.a      = para.out$coef.a,
        coef.b      = para.out$coef.b,
        Theta       = para.out$Theta,
        X1.cor.list = para.out$X1.cor
      )
    }
    
    if (print.detail) {
      message("Worker [", Sys.getpid(), "] finished community ", i)
    }
    
    list(community = i, result = local_para)
  })
      
      # Update corrected data (X.cor.1) for this non-reference batch
      # We update for each valid group; note that the order in Xb1.Batk.COMi.glist now corresponds 
      # to the subset of groups with data (which we can retrieve from valid_groups_nonref)
  for (res in comm_res) {
    if (is.null(res)) next
    for (k in seq_along(res$result)) {
      entry <- res$result[[k]]
      if (is.null(entry)) next
      
      metID  <- entry$metID
      v1     <- entry$valid1
      ga     <- which(v1)           # positions in group.levels
      X1list <- entry$X1.cor.list
      
      ## 1) fill X.cor.1
      for (j in seq_along(ga)) {
        gpos <- ga[j]
        rows <- which(batch == batch.levels[k] & group == group.levels[gpos])
        X.cor.1[rows, metID] <- X1list[[j]]
      }
      
      ## 2) record parameters
      Xcor.para[[k]]$coef.a[metID] <- entry$coef.a
      Xcor.para[[k]]$coef.b[metID] <- entry$coef.b
      for (j in seq_along(ga)) {
        Xcor.para[[k]]$Theta[[ ga[j] ]][ metID, metID ] <- entry$Theta[[j]]
      }
    }
  }
  
  X.cor <- matrix(0, nrow = nrow(X.nodel), ncol = p)
  for (k in seq_len(batch.num)) {
    lbl <- batch.levels[k]
    rows <- which(batch.old == lbl)
    if (length(rows)==0) next
    A <- diag(Xcor.para[[k]]$coef.a)
    B <- matrix(Xcor.para[[k]]$coef.b, nrow = length(rows), ncol = p, byrow = TRUE)
    X.cor[rows, ] <- X.nodel[ rows, ] %*% A + B
  }
  
  ## --- build X.cor.withQC if needed ---
  X.cor.withQC <- NULL
  if (containQC) {
    X.cor.withQC <- matrix(0, nrow = nrow(X.init), ncol = p)
    ## ref‐batch QCs and non‐QC rows
    refq <- which(batch.init == ref.batch_char & group.init=="QC")
    if (length(refq)) X.cor.withQC[refq, ] <- X.init[refq, ]
    nonQC <- which(group.init != "QC")
    X.cor.withQC[nonQC, ] <- X.cor
  }
  
    
#    if (print.detail) message("Finished correction of community ", i)
#  }
  
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
