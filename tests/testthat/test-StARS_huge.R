library(testthat)
library(CordBat)  

test_that("StARS_huge returns a named numeric vector of length 2", {
  skip_if_not_installed("huge")
  set.seed(123)
  
  # small synthetic matrix
  X <- matrix(rnorm(200), nrow = 20, ncol = 10)
  res <- CordBat:::StARS_huge(
    X       = X,
    b       = round(0.7 * nrow(X)),
    M       = 5,
    lambda.grid = seq(0.9,0.1,-0.1),
    seed = NULL,
    beta    = 0.05,
    print.detail = FALSE
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
    b       = round(0.7 * nrow(X)),
    M       = 5,
    lambda.grid = seq(0.9,0.1,-0.1),
    seed = NULL,
    beta    = 0.05,
    print.detail = FALSE
  )
  
  # selected rho should be positive
  expect_gt(res["Sel.rho"], 0)
  # instability must lie in [0,1]
  expect_true(res["D_var"] >= 0 && res["D_var"] <= 1)
})

test_that("print.detail = TRUE emits a selection message", {
  skip_if_not_installed("huge")
  set.seed(789)
  
  X <- matrix(rnorm(200), nrow = 20, ncol = 10)
  expect_message(
    CordBat:::StARS_huge(
      X       = X,
      b       = round(0.7 * nrow(X)),
      M       = 5,
      lambda.grid = seq(0.9,0.1,-0.1),
      seed = NULL,
      beta    = 0.05,
      print.detail = TRUE
    ),
    "StARS selected rho"
  )
})