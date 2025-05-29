library(testthat)
library(CordBat)

test_that("select_rho_cv_bic returns correct structure and values for CV case", {
  skip_if_not_installed("CVglasso")
  skip_if_not_installed("huge")
  
  # create a data set large enough so selfoldforCV() picks fold > 1
  X <- matrix(rnorm(50 * 5), nrow = 50, ncol = 5)
  rhos <- c(0.1, 0.5, 0.9)
  
  res <- CordBat:::select_rho_cv_bic(X, rhos,
                           print.detail = FALSE,
                           seed         = 123)
  
  # named numeric vector of length 2
  expect_named(res, c("rho.sel", "MinCVerr"))
  expect_length(res, 2)
  expect_type(res["rho.sel"],   "double")
  expect_type(res["MinCVerr"],  "double")
  
  # rho.sel must be one of the inputs
  expect_true(res["rho.sel"] %in% rhos)
})

test_that("select_rho_cv_bic returns correct structure and values for BIC fallback", {
  skip_if_not_installed("CVglasso")
  skip_if_not_installed("huge")
  
  # small data so selfoldforCV() returns fold == 1
  X <- matrix(rnorm(5 * 5), nrow = 5, ncol = 5)
  rhos <- seq(0.1, 0.9, by = 0.1)
  
  res <- CordBat:::select_rho_cv_bic(X, rhos,
                           print.detail = FALSE,
                           seed         = 123)
  
  expect_named(res, c("rho.sel", "MinCVerr"))
  expect_length(res, 2)
  expect_type(res["rho.sel"],   "double")
  expect_type(res["MinCVerr"],  "double")
  
  expect_true(res["rho.sel"] %in% rhos)
})

test_that("print.detail = TRUE emits the correct message", {
  skip_if_not_installed("CVglasso")
  skip_if_not_installed("huge")
  
  X <- matrix(rnorm(50 * 5), nrow = 50, ncol = 5)
  rhos <- c(0.1, 0.5, 0.9)
  
  expect_message(
    CordBat:::select_rho_cv_bic(X, rhos, print.detail = TRUE, seed = 123),
    "CV \\+ BIC selects rho ="
  )
})
