library(testthat)
library(CordBat)

test_that("graphicalLasso delegates exactly to glassoFast", {
  set.seed(42)
  X   <- matrix(rnorm(200), nrow = 40, ncol = 5)
  rho <- 0.15
  
  # mirror what your function does internally
  Xs <- scale(X, center = TRUE, scale = TRUE)
  S  <- cov(Xs)
  
  rho.mat <- matrix(rho, nrow(S), ncol(S))
  diag(rho.mat) <- 0
  # 1) call your new wrapper
  out_wrap <- CordBat:::graphicalLasso(X, rho, print.detail = FALSE)
  
  # 2) call glassoFast directly
  out_fast <- glassoFast::glassoFast(S, rho.mat, thr = 1e-5)
  
  # they should match to floating‐point tolerance
  expect_equal(out_wrap$W,     out_fast$w,  tolerance = 1e-8)
  expect_equal(out_wrap$Theta, out_fast$wi, tolerance = 1e-8)
})

test_that("rho = 0 returns raw covariance and exact inverse (within 1e-5)", {
  set.seed(42)
  X <- matrix(rnorm(300), nrow = 60, ncol = 5)
  res0    <- CordBat:::graphicalLasso(X, rho = 0, print.detail = FALSE)
  S0      <- cov(scale(X, TRUE, TRUE))
  Theta0  <- solve(S0)
  
  dW      <- abs(res0$W     - S0)
  dTheta  <- abs(res0$Theta - Theta0)
  tol     <- 1e-5
  
  expect_true(max(dW)     < tol, info = paste("max |W - S0|     =", max(dW)))
  expect_true(max(dTheta) < tol, info = paste("max |Θ - inv(S0)| =", max(dTheta)))
})