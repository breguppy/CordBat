# tests/testthat/test-CordBat.R

library(testthat)
library(CordBat)

# If CordBat is part of a package, replace with the correct library() call.
# For example, if CordBat is exported by a package named "MyPackage", uncomment the line below:
# library(MyPackage)

test_that("CordBat outputs have correct structure with skip.impute = TRUE", {
  set.seed(123)
  X <- matrix(rnorm(20), nrow = 5, ncol = 4)
  batch <- rep(c("A", "B"), length.out = 5)
  
  res <- CordBat:::CordBat(
    X          = X,
    batch      = batch,
    ref.batch  = "A",
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
  expect_equal(res$batch.new, as.character(batch))
  
  # 4. X.delout should be identical to the original X (no outliers removed)
  expect_true(is.matrix(res$X.delout))
  expect_equal(dim(res$X.delout), dim(X))
  expect_equal(res$X.delout, X)
  
  # 5. X.cor should be a matrix with the same number of columns as X
  expect_true(is.matrix(res$X.cor))
  expect_equal(ncol(res$X.cor), ncol(X))
  
  # 6. For rows in the reference batch ("A"), X.cor must equal the original X
  ref_idx <- which(batch == "A")
  expect_equal(res$X.cor[ref_idx, ], X[ref_idx, ])
  
  # 7. X.cor.1 (correction on X.delout) for the reference batch should also match X
  expect_true(is.matrix(res$X.cor.1))
  expect_equal(ncol(res$X.cor.1), ncol(X))
  expect_equal(res$X.cor.1[ref_idx, ], X[ref_idx, ])
  
  # 8. Because there were no QC samples in this test, X.cor.withQC should be NULL
  expect_null(res$X.cor.withQC)
})

test_that("CordBat handles QC samples when grouping = TRUE", {
  set.seed(456)
  # Construct a dataset with QC samples such that the reference batch has >=2 non-QC samples
  X <- matrix(rnorm(24), nrow = 6, ncol = 4)
  batch <- rep(c("A", "B"), each = 3)
  # Assign groups: for batch "A": "1", "QC", "2"; for batch "B": "QC", "3", "QC"
  group <- c("1", "QC", "2", "QC", "3", "QC")
  
  res <- CordBat:::CordBat(
    X         = X,
    batch     = batch,
    group     = group,
    grouping  = TRUE,
    ref.batch = "A",
    skip.impute = TRUE,
    print.detail = FALSE
  )
  
  # 1. Because there are QC samples, X.cor.withQC must be a matrix of the same dimensions as X
  expect_true(is.matrix(res$X.cor.withQC))
  expect_equal(dim(res$X.cor.withQC), dim(X))
  
  # 2. For QC samples in the reference batch ("A"), X.cor.withQC should equal the original X
  ref_qc_idx <- which(batch == "A" & group == "QC")
  expect_equal(res$X.cor.withQC[ref_qc_idx, ], X[ref_qc_idx, ])
  
  # 3. For non-QC samples in the reference batch, X.cor should equal the original X
  nonQC_idx <- which(group != "QC")
  ref_nonqc_orig <- which(batch == "A" & group != "QC")
  ref_nonqc_positions <- match(ref_nonqc_orig, nonQC_idx)
  
  expect_equal(
    res$X.cor[ref_nonqc_positions, ],
    X[ref_nonqc_orig, ]
  )
  
  # 4. Because skip.impute = TRUE, no samples are removed
  expect_equal(res$delsampIdx, integer(0))
})
