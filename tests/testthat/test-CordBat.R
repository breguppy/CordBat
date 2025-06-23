library(testthat)
library(CordBat)

sim_data <- function(n_batch = 2, per_batch = 6, p = 4, add_qc = FALSE) {
  set.seed(333)
  n  <- n_batch * per_batch + if (add_qc) 2 else 0
  X  <- matrix(rnorm(n * p), n, p)
  batch <- rep(seq_len(n_batch), each = per_batch)
  if (add_qc) {
    group <- c(rep(1, n - 2), "QC", "QC")
    batch[(n - 1):n] <- 1                      # QC in ref batch
  } else {
    group <- rep(1, n)
  }
  list(X = X, batch = batch, group = group)
}

test_that("CordBat outputs have correct structure with skip.impute = TRUE", {
  dat <- sim_data(n_batch = 3, per_batch = 6, p = 10, add_qc = FALSE)
  
  res <- CordBat(
    X          = dat$X,
    batch      = dat$batch,
    ref.batch  = 1,
    skip.impute = TRUE,
    print.detail = FALSE
  )
  
  # 1. The result should be a list with exactly these names
  expected_names <- c(
    "batch.level", "delsampIdx", "batch.new", "group.new",
    "X.delout", "X.cor", "X.cor.1", "X.cor.withQC", "Xcor.para"
  )
  expect_type(res, "list")
  expect_named(res, expected_names)
  
  # 2. Since skip.impute = TRUE, no samples are deleted
  expect_equal(res$delsampIdx, integer(0))
  
  # 3. batch.new should be the same as the input, coerced to character
  expect_equal(res$batch.new, as.character(dat$batch))
  
  # 4. X.delout should be identical to the original X (no outliers removed)
  expect_true(is.matrix(res$X.delout))
  expect_equal(dim(res$X.delout), dim(dat$X))
  expect_equal(res$X.delout, dat$X)
  
  # 5. X.cor should be a matrix with the same number of columns as X
  expect_true(is.matrix(res$X.cor))
  expect_equal(ncol(res$X.cor), ncol(dat$X))
  
  # 6. For rows in the reference batch, X.cor must equal the original X
  ref_idx <- which(dat$batch == 1)
  expect_equal(res$X.cor[ref_idx, ], dat$X[ref_idx, ])
  
  # 7. X.cor.1 (correction on X.delout) for the reference batch should also match X
  expect_true(is.matrix(res$X.cor.1))
  expect_equal(ncol(res$X.cor.1), ncol(dat$X))
  expect_equal(res$X.cor.1[ref_idx, ], dat$X[ref_idx, ])
  
  # 8. Because there were no QC samples in this test, X.cor.withQC should be NULL
  expect_null(res$X.cor.withQC)
  
})

test_that("CordBat handles QC samples when grouping = TRUE", {
  dat <- sim_data(add_qc = TRUE)
  res <- CordBat(dat$X, dat$batch, dat$group,
                 grouping     = TRUE,
                 ref.batch    = 1,
                 print.detail = FALSE,
                 skip.impute  = TRUE)
  
  # 1. Because there are QC samples, X.cor.withQC must be a matrix of the same dimensions as X
  expect_true(is.matrix(res$X.cor.withQC))
  expect_equal(dim(res$X.cor.withQC), dim(dat$X))
  
  # 2. For QC samples in the reference batch, X.cor.withQC should equal the original X
  ref_qc_idx <- which(batch == 1 & group == "QC")
  expect_equal(res$X.cor.withQC[ref_qc_idx, ], dat$X[ref_qc_idx, ])
  
  # 3. For non-QC samples in the reference batch, X.cor should equal the original X
  nonQC_idx <- which(dat$group != "QC")
  ref_nonqc_orig <- which(dat$batch == 1 & dat$group != "QC")
  ref_nonqc_positions <- match(ref_nonqc_orig, nonQC_idx)
  
  expect_equal(
    res$X.cor[ref_nonqc_positions, ],
    dat$X[ref_nonqc_orig, ]
  )
  
  # 4. Because skip.impute = TRUE, no samples are removed
  expect_equal(res$delsampIdx, integer(0))
})

test_that("CordBat removes an extreme outlier from X.cor.1 when skip.impute = FALSE", {
  set.seed(42)
  # Create a 10×4 matrix with one extreme outlier in the last row
  X_normal  <- matrix(rnorm(10 * 4, mean = 0, sd = 1), nrow = 10, ncol = 4)
  X_outlier <- matrix(1e5, nrow = 1, ncol = 4)
  X <- rbind(X_normal, X_outlier)
  
  batch <- rep("A", 11)
  
  res_outlier <- DelOutlier(X)
  n_removed_outlier <- length(res_outlier$delsampIdx)
  
  expect_gt(n_removed_outlier, 0)
  
  res <- CordBat(
    X           = X,
    batch       = batch,
    ref.batch   = "A",
    skip.impute = FALSE,
    print.detail = FALSE
  )
  
  # The extreme outlier is the 11th sample, so delsampIdx should be exactly 11
  expect_length(res$delsampIdx, 1)
  expect_equal(res$delsampIdx, 11)
  
  # X.cor always has the same number of rows as the original X (11 rows)
  expect_equal(nrow(res$X.cor), 11)
  
  # X.cor.1 should have one fewer row because the extreme outlier was removed
  expect_lt(nrow(res$X.cor.1), nrow(res$X.cor))
  
  # X.cor.1 should have one fewer row than X
  expect_true(is.matrix(res$X.cor.1))
  expect_equal(nrow(res$X.cor.1), nrow(X) - 1)
  
  # None of the rows in X.cor.1 should match the extreme vector c(1e6,1e6,1e6,1e6)
  extreme_vec <- as.numeric(X_outlier)
  found_matching <- any(apply(res$X.cor.1, 1, function(row) all(row == extreme_vec)))
  expect_false(found_matching)
})

test_that("CordBat: print.detail produces messages", {
  dat <- sim_data()
  #expect_silent(
  #  CordBat(dat$X, dat$batch, ref.batch = 1,
  #          print.detail = FALSE, skip.impute = TRUE)
  #)
  
  expect_message(
    CordBat(dat$X, dat$batch, ref.batch = 1,
            print.detail = TRUE, skip.impute = TRUE),
    regexp = ".", all = FALSE
  )
})


test_that("CordBat errors if ref.batch is not present", {
  dat <- sim_data()
  expect_error(
    CordBat(dat$X, dat$batch, ref.batch = "nonexistent",
            print.detail = FALSE, skip.impute = TRUE),
    regexp = "invalid"
  )
})

test_that("CordBat correctly saves correction parameters", {
  dat <- sim_data(n_batch = 3, per_batch = 10, p = 15, add_qc = FALSE)
  
  res <- CordBat(
    X          = dat$X,
    batch      = dat$batch,
    ref.batch  = 1,
    skip.impute = TRUE,
    print.detail = FALSE
  )
  
  expect_type(res$Xcor.para, "list")
  
  # For rows in the reference batch, X.cor must equal the original X 
  ref_idx <- which(dat$batch == 1)
  expect_equal(res$X.cor[ref_idx, ], dat$X[ref_idx, ])
  # coef.a should be all 1 and coef.b should be all 0
  expect_equal(res$Xcor.para[[1]]$coef.a, rep(1,15))
  expect_equal(res$Xcor.para[[1]]$coef.b, rep(0,15))
  
  # corrected batch 2 and 3 should equal A * X + b
  batch2_idx <- which(dat$batch == 2)
  batch2_Amat <- diag(res$Xcor.para[[2]]$coef.a)
  batch2_bmat <- matrix(rep(res$Xcor.para[[2]]$coef.b, each = 10), nrow = 10)
  batch2_cor <- dat$X[batch2_idx, ] %*% batch2_Amat + batch2_bmat
  expect_equal(res$X.cor[batch2_idx, ], batch2_cor)
  
  batch3_idx <- which(dat$batch == 3)
  batch3_Amat <- diag(res$Xcor.para[[3]]$coef.a)
  batch3_bmat <- matrix(rep(res$Xcor.para[[3]]$coef.b, each = 10), nrow = 10)
  batch3_cor <- dat$X[batch3_idx, ] %*% batch3_Amat + batch3_bmat
  expect_equal(res$X.cor[batch3_idx, ], batch3_cor)
})

test_that("CordBat: non-QC rows in X.cor.withQC match X.cor", {
  # 1) Simulate a small toy dataset
  set.seed(42)
  X <- matrix(rnorm(8*3), nrow = 8, ncol = 3,
              dimnames = list(NULL, paste0("M", 1:3)))
  
  # 2 batches of 4 each; in each batch mark the 4th row as QC
  batch <- rep(1:2, each = 4)
  group <- c("A", "A", "A", "QC", "B", "B", "B", "QC")
  
  # 2) Run CordBat, skipping outlier imputation for speed
  res <- CordBat(
    X           = X,
    batch       = batch,
    group       = group,
    ref.batch   = 1,
    skip.impute = TRUE,
    print.detail= FALSE
  )
  
  # 3) Identify non-QC indices
  nonQC_orig_idx <- which(group != "QC")  
  
  # 4) Assert exact equality
  expect_equal(
    res$X.cor.withQC[nonQC_orig_idx, , drop = FALSE],
    res$X.cor
  )
})
