# tests/testthat/test-BEgLasso.R

library(testthat)
library(CordBat)

test_that("BEgLasso returns a list with the correct components and dimensions", {
  set.seed(42)
  # two groups, each with a 2×2 matrix
  X0.glist <- list(matrix(rnorm(4), 2, 2), matrix(rnorm(4), 2, 2))
  X1.glist <- list(matrix(rnorm(4), 2, 2), matrix(rnorm(4), 2, 2))
  
  res <- CordBat:::BEgLasso(
    X0.glist, X1.glist,
    penal.rho   = 0.1,
    penal.ksi   = 0.1,
    penal.gamma = 0.1,
    eps         = 1e-4,
    print.detail = FALSE
  )
  
  # 1) Top‐level structure
  expect_type(res, "list")
  expect_named(res, c("Theta", "X1.cor", "coef.a", "coef.b"))
  
  # 2) Theta: list of G precision matrices, each p×p
  expect_length(res$Theta, length(X0.glist))
  for (theta in res$Theta) {
    expect_true(is.matrix(theta))
    expect_equal(dim(theta), c(2, 2))
  }
  
  # 3) X1.cor: list of corrected matrices, same dims as inputs
  expect_length(res$X1.cor, length(X1.glist))
  for (mat in res$X1.cor) {
    expect_true(is.matrix(mat))
    expect_equal(dim(mat), dim(X1.glist[[1]]))
  }
  
  # 4) coef.a and coef.b: numeric vectors of length p, coef.a ≥ 0
  expect_numeric(res$coef.a, any.missing = FALSE, len = 2)
  expect_true(all(res$coef.a >= 0))
  expect_numeric(res$coef.b, any.missing = FALSE, len = 2)
})

test_that("With identical X0 and X1 and zero penalties, no correction occurs", {
  set.seed(42)
  X <- matrix(rnorm(4), 2, 2)
  X0.glist <- list(X, X)
  X1.glist <- list(X, X)
  
  res <- CordBat:::BEgLasso(
    X0.glist, X1.glist,
    penal.rho   = 0,
    penal.ksi   = 0,
    penal.gamma = 0,
    eps         = 1e-6,
    print.detail = FALSE
  )
  
  # For zero penalization and identical data, corrected = original
  expect_equal(res$X1.cor, X1.glist)
})