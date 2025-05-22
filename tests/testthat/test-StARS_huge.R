library(testthat)
library(CordBat)  

test_that("StARS_huge returns a named numeric vector of length 2", {
  skip_if_not_installed("huge")
  set.seed(123)
  
  # small synthetic matrix
  X <- matrix(rnorm(200), nrow = 20, ncol = 10)
  res <- CordBat:::StARS_huge(
    X       = X,
    b       = ceiling(0.7 * nrow(X)),
    M       = 5,
    nlambda = 10,
    beta    = 0.1,
    verbose = FALSE
  )
  
  expect_type(res, "double")
  expect_length(res, 2)
  expect_named(res, c("Sel.rho", "D_var"))
})

test_that("Returned values are in plausible ranges", {
  skip_if_not_installed("huge")
  set.seed(456)
  
  X <- matrix(rnorm(300), nrow = 30, ncol = 10)
  res <- CordBat:::StARS_huge(
    X       = X,
    b       = ceiling(0.7 * nrow(X)),
    M       = 5,
    nlambda = 10,
    beta    = 0.2,
    verbose = FALSE
  )
  
  # selected rho should be positive
  expect_gt(res["Sel.rho"], 0)
  # instability must lie in [0,1]
  expect_true(res["D_var"] >= 0 && res["D_var"] <= 1)
})

test_that("Verbose = TRUE emits a selection message", {
  skip_if_not_installed("huge")
  set.seed(789)
  
  X <- matrix(rnorm(200), nrow = 20, ncol = 10)
  expect_message(
    CordBat:::StARS_huge(
      X       = X,
      b       = ceiling(0.7 * nrow(X)),
      M       = 3,
      nlambda = 5,
      beta    = 0.3,
      verbose = TRUE
    ),
    "StARS selected rho"
  )
})