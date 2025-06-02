library(testthat)
library(mixOmics)  # Needed for pca
library(CordBat)  # Adjust to correct path if needed

test_that("DelOutlier detects and removes PCA outliers", {
  set.seed(123)
  X <- matrix(rnorm(300), nrow = 30, ncol = 10) 
  
  # Run without outliers
  res_clean <- CordBat:::DelOutlier(X)
  n_removed_clean <- length(res_clean$delsampIdx)
  
  # Add a clear outlier to row 1
  X_with_outlier <- X
  X_with_outlier[1, ] <- X_with_outlier[1, ] + 20
  
  res_outlier <- DelOutlier(X_with_outlier)
  n_removed_outlier <- length(res_outlier$delsampIdx)
  
  # Ensure outlier was removed
  expect_gt(n_removed_outlier, n_removed_clean)
  expect_true(1 %in% res_outlier$delsampIdx)
  
  # Ensure the output structure is intact
  expect_named(res_outlier, c("delsampIdx", "X.out"))
  expect_equal(ncol(res_outlier$X.out), ncol(X))
  expect_equal(nrow(res_outlier$X.out), nrow(X) - n_removed_outlier)
})
