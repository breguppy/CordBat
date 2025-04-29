library(testthat)
library(CordBat)

test_that("CDfgL converges on simple input", {
  V <- diag(3)  # Identity matrix
  beta_i <- c(0.5, -0.3, 0.1)
  u <- c(1, 1, 1)
  rho <- 0.1
  
  result <- CordBat:::CDfgL(V, beta_i, u, rho, maxIter = 100, print.detail = FALSE)
  
  expect_type(result, "double")
  expect_length(result, 3)
  expect_true(all(is.finite(result)))
})

test_that("CDfgL returns zero for zero input", {
  V <- diag(3)
  beta_i <- c(0, 0, 0)
  u <- c(0, 0, 0)
  rho <- 0.5
  
  result <- CordBat:::CDfgL(V, beta_i, u, rho, maxIter = 50, print.detail = FALSE)
  
  expect_equal(result, c(0, 0, 0))
})

test_that("CDfgL hits maximum iterations when convergence is impossible", {
  set.seed(123)
  V <- diag(c(1, 2, 3))  # Diagonal matrix but not identity
  beta_i <- runif(3, -10, 10)  # Start far from zero
  u <- runif(3, -10, 10)
  rho <- 1e-10  # Tiny regularization -> updates will be very small
  
  expect_message(
    CordBat:::CDfgL(V, beta_i, u, rho, maxIter = 1, print.detail = TRUE),
    "CDfgL reached max iteration"
  )
})