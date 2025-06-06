library(testthat)
library(CordBat)

# ------------------------------------------------------------------
# Helper to simulate quick toy matrices ----------------------------
# ------------------------------------------------------------------
simulate_X <- function(n = 20, p = 7) {
  matrix(rnorm(n * p), n, p)
}

test_that("graphicalLasso returns expected top‑level structure", {
  set.seed(11)
  X <- simulate_X(18, 6)
  
  res <- graphicalLasso(X, rho = 0.1, print.detail = FALSE)
  
  expect_type(res, "list")
  expect_named(res, c("Theta", "W"), ignore.order = TRUE)
  
  # Dimensions & symmetry ------------------------------------------
  expect_equal(dim(res$W), dim(res$Theta))
  expect_equal(res$Theta, t(res$Theta), tolerance = 1e-10)
  expect_equal(res$W, t(res$W), tolerance = 1e-10)
})

test_that("graphicalLasso Precision and covariance matrices are positive‑definite", {
  set.seed(202)
  X <- simulate_X(25, 8)
  out <- graphicalLasso(X, rho = 0.15, print.detail = FALSE)
  
  eigW     <- eigen(out$W, symmetric = TRUE, only.values = TRUE)$values
  eigTheta <- eigen(out$Theta, symmetric = TRUE, only.values = TRUE)$values
  expect_true(all(eigW > 0))
  expect_true(all(eigTheta > 0))
})

test_that("graphicalLasso print.detail toggles console output / messages", {
  X <- simulate_X(12, 5)
  
  expect_silent(
    graphicalLasso(X, rho = 0.05, print.detail = FALSE)
  )
  
})

test_that("graphicalLasso delegates exactly to glassoFast", {
  set.seed(42)
  X   <- simulate_X(40, 5)
  rho <- 0.15
  
  # mirror what your function does internally
  Xs <- scale(X, center = TRUE, scale = TRUE)
  S  <- cov(Xs)
  
  rho.mat <- matrix(rho, nrow(S), ncol(S))
  #diag(rho.mat) <- 0
  # 1) call your new wrapper
  out_wrap <- graphicalLasso(X, rho, print.detail = FALSE)
  
  # 2) call glassoFast directly
  out_fast <- glassoFast::glassoFast(S, rho.mat, thr = 1e-5)
  
  # they should match to floating‐point tolerance
  expect_equal(out_wrap$W,     out_fast$w,  tolerance = 1e-8)
  expect_equal(out_wrap$Theta, out_fast$wi, tolerance = 1e-8)
})

test_that("graphicalLasso rho = 0 returns raw covariance and exact inverse (within 1e-5)", {
  set.seed(42)
  X <- simulate_X(60, 5)
  res0    <- graphicalLasso(X, rho = 0, print.detail = FALSE)
  S0      <- cov(scale(X, TRUE, TRUE))
  Theta0  <- solve(S0)
  
  dW      <- abs(res0$W - S0)
  dTheta  <- abs(res0$Theta - Theta0)
  tol     <- 1e-5
  
  expect_true(max(dW) < tol, info = paste("max |W - S0| =", max(dW)))
  expect_true(max(dTheta) < tol, info = paste("max |Theta - inv(S0)| =", max(dTheta)))
})

test_that("graphicalLasso produces same output as old_graphicalLasso", {
  set.seed(456)
  X <- simulate_X(50, 10)
  
  new_res <- graphicalLasso(X, rho = 0.9, print.detail = FALSE)
  old_res <- old_graphicalLasso(X, rho = 0.9, print.detail = FALSE)
  
  expect_equal(new_res$W, old_res$W,  tolerance = 1e-8)
  expect_equal(new_res$Theta, old_res$Theta, tolerance = 1e-8)
  
})

