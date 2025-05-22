library(testthat)
library(CordBat)

test_that("selrho.useCVBIC returns a numeric vector of length 2", {
  skip_if_not_installed("huge")
  set.seed(123)
  
  X <- matrix(rnorm(200), nrow = 20, ncol = 10)
  res <- CordBat:::selrho.useCVBIC(
    X            = X,
    print.detail = FALSE,
    nlambda      = 10,
    gamma        = 0
  )
  
  expect_type(res, "double")
  expect_length(res, 2)
})

test_that("returned values are in plausible ranges", {
  skip_if_not_installed("huge")
  set.seed(456)
  
  X <- matrix(rnorm(300), nrow = 30, ncol = 10)
  res <- CordBat:::selrho.useCVBIC(
    X            = X,
    print.detail = FALSE,
    nlambda      = 10,
    gamma        = 0.5
  )
  
  # first element is the selected rho > 0
  expect_gt(res[1], 0)
  # second element is the EBIC score (finite, non-NA)
  expect_true(is.finite(res[2]))
})

test_that("verbose = TRUE emits an EBIC selection message", {
  skip_if_not_installed("huge")
  set.seed(789)
  
  X <- matrix(rnorm(200), nrow = 20, ncol = 10)
  expect_message(
    CordBat:::selrho.useCVBIC(
      X            = X,
      print.detail = TRUE,
      nlambda      = 5,
      gamma        = 0
    ),
    "EBIC: select rho"
  )
})
