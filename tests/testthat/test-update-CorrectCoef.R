library(testthat)
library(CordBat)

test_that("update.CorrectCoef returns the expected datatypes and format", {
  
  # Simple mock data
  set.seed(123)
  
  # 2 groups, each with 3 samples, 2 features
  X0.glist <- list(
    matrix(rnorm(6), nrow = 3, ncol = 2),
    matrix(rnorm(6), nrow = 3, ncol = 2)
  )
  
  X1.glist <- list(
    matrix(rnorm(6), nrow = 3, ncol = 2),
    matrix(rnorm(6), nrow = 3, ncol = 2)
  )
  
  # Theta matrices (precision-like matrices), size 2x2 for each group
  Theta.list <- list(
    diag(2), 
    diag(2)
  )
  
  # Initial a and b coefficients
  a.i <- rep(1, 2)
  b.i <- rep(0, 2)
  
  penal.ksi <- 0.1
  penal.gamma <- 0.1
  print.detail <- FALSE
  
  # Call the function
  result <- CordBat:::update.CorrectCoef(
    X0.glist = X0.glist,
    X1.glist = X1.glist,
    Theta.list = Theta.list,
    a.i = a.i,
    b.i = b.i,
    penal.ksi = penal.ksi,
    penal.gamma = penal.gamma,
    print.detail = print.detail
  )
  
  # Check output structure
  expect_type(result, "list")
  expect_named(result, c("coef.a", "coef.b"))
  
  # Check length of output vectors
  expect_length(result$coef.a, 2)
  expect_length(result$coef.b, 2)
  
  # Optional: Check numerical values (basic sanity checks)
  expect_true(all(is.finite(result$coef.a)))
  expect_true(all(is.finite(result$coef.b)))
})

test_that("update.CorrectCoef updates coefficients correctly with fixed input", {
  
  # Fixed mock data
  X0.glist <- list(
    matrix(c(1, 2, 3, 4), nrow = 2, ncol = 2, byrow = TRUE)
  )
  
  X1.glist <- list(
    matrix(c(5, 6, 7, 8), nrow = 2, ncol = 2, byrow = TRUE)
  )
  
  Theta.list <- list(
    diag(2)
  )
  
  a.i <- c(1, 1)
  b.i <- c(0, 0)
  
  penal.ksi <- 0.1
  penal.gamma <- 0.1
  print.detail <- FALSE
  
  result <- CordBat:::update.CorrectCoef(
    X0.glist = X0.glist,
    X1.glist = X1.glist,
    Theta.list = Theta.list,
    a.i = a.i,
    b.i = b.i,
    penal.ksi = penal.ksi,
    penal.gamma = penal.gamma,
    print.detail = print.detail
  )
  
  # Now match the actual observed output
  expected_coef_a <- c(0.67, 0.71)
  expected_coef_b <- c(-0.77, -0.76)
  
  expect_equal(result$coef.a, expected_coef_a, tolerance = 1e-2)
  expect_equal(result$coef.b, expected_coef_b, tolerance = 1e-2)
  
  # Also check structure
  expect_named(result, c("coef.a", "coef.b"))
  expect_length(result$coef.a, 2)
  expect_length(result$coef.b, 2)
})

test_that("update.CorrectCoef updates coefficients reasonably", {
  
  # Fixed mock data
  X0.glist <- list(
    matrix(c(1, 2, 3, 4), nrow = 2, ncol = 2, byrow = TRUE)
  )
  
  X1.glist <- list(
    matrix(c(5, 6, 7, 8), nrow = 2, ncol = 2, byrow = TRUE)
  )
  
  Theta.list <- list(
    diag(2)
  )
  
  a.i <- c(1, 1)
  b.i <- c(0, 0)
  
  penal.ksi <- 0.1
  penal.gamma <- 0.1
  print.detail <- FALSE
  
  result <- CordBat:::update.CorrectCoef(
    X0.glist = X0.glist,
    X1.glist = X1.glist,
    Theta.list = Theta.list,
    a.i = a.i,
    b.i = b.i,
    penal.ksi = penal.ksi,
    penal.gamma = penal.gamma,
    print.detail = print.detail
  )
  
  # --- Behavioral Tests ---
  
  # Structure checks
  expect_named(result, c("coef.a", "coef.b"))
  expect_length(result$coef.a, 2)
  expect_length(result$coef.b, 2)
  
  # Value behavior checks
  
  # a coefficients should remain positive
  expect_true(all(result$coef.a > 0))
  
  # a coefficients should shrink or stay near 1 (no big explosion)
  expect_true(all(result$coef.a <= 1))
  
  # b coefficients should be finite
  expect_true(all(is.finite(result$coef.b)))
  
  # b coefficients might become slightly negative
  expect_true(all(result$coef.b < 0))
  
  # coefficients should not be extreme
  expect_true(all(abs(result$coef.a) < 10))
  expect_true(all(abs(result$coef.b) < 10))
})

test_that("update.CorrectCoef handles Sigma_g[j] == 0 safely", {
  warnings <- NULL
  
  result <- withCallingHandlers(
    CordBat:::update.CorrectCoef(
      X0.glist = list(matrix(c(1, 1, 2, 1), nrow = 2, byrow = TRUE)),
      X1.glist = list(matrix(c(3, 1, 4, 1), nrow = 2, byrow = TRUE)),
      Theta.list = list(diag(2)),
      a.i = c(1, 1),
      b.i = c(0, 0),
      penal.ksi = 0.1,
      penal.gamma = 0.1,
      print.detail = TRUE
    ),
    warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  
  # Check that a warning about Sigma_g was captured
  expect_true(any(grepl("Sigma_g\\[\\d+\\] was zero or NA", warnings)))
  
  # Check output structure and values
  expect_named(result, c("coef.a", "coef.b"))
  expect_length(result$coef.a, 2)
  expect_length(result$coef.b, 2)
  expect_true(all(is.finite(result$coef.a)))
  expect_true(all(is.finite(result$coef.b)))
})

test_that("update.CorrectCoef handles multiple groups with Sigma_g[j] == 0", {
  warnings <- NULL
  
  X0.glist <- list(
    matrix(c(1, 1, 2, 1), nrow = 2, byrow = TRUE),  # col2 is constant
    matrix(c(2, 2, 3, 2), nrow = 2, byrow = TRUE)   # col2 is constant
  )
  
  X1.glist <- list(
    matrix(c(3, 1, 4, 1), nrow = 2, byrow = TRUE),  # col2 is constant
    matrix(c(5, 2, 6, 2), nrow = 2, byrow = TRUE)   # col2 is constant
  )
  
  Theta.list <- list(diag(2), diag(2))
  
  result <- withCallingHandlers(
    CordBat:::update.CorrectCoef(
      X0.glist = X0.glist,
      X1.glist = X1.glist,
      Theta.list = Theta.list,
      a.i = c(1, 1),
      b.i = c(0, 0),
      penal.ksi = 0.1,
      penal.gamma = 0.1,
      print.detail = TRUE
    ),
    warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  
  expect_true(length(warnings) >= 2)
  expect_true(all(grepl("Sigma_g\\[\\d+\\] was zero or NA", warnings)))
  expect_named(result, c("coef.a", "coef.b"))
  expect_true(all(is.finite(result$coef.a)))
  expect_true(all(is.finite(result$coef.b)))
})

test_that("update.CorrectCoef handles Sigma_g[j] == NA safely", {
  warnings <- NULL
  
  # All NA in second column → scale() will return NA for scale factor
  X0.glist <- list(
    matrix(c(1, NA, 2, NA), nrow = 2, byrow = TRUE)
  )
  
  X1.glist <- list(
    matrix(c(3, NA, 4, NA), nrow = 2, byrow = TRUE)
  )
  
  Theta.list <- list(diag(2))
  
  result <- withCallingHandlers(
    CordBat:::update.CorrectCoef(
      X0.glist = X0.glist,
      X1.glist = X1.glist,
      Theta.list = Theta.list,
      a.i = c(1, 1),
      b.i = c(0, 0),
      penal.ksi = 0.1,
      penal.gamma = 0.1,
      print.detail = TRUE
    ),
    warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  
  expect_true(any(grepl("Sigma_g\\[\\d+\\] was zero or NA", warnings)))
  expect_true(all(is.finite(result$coef.a)))
  expect_true(all(is.finite(result$coef.b)))
})
test_that("C++ core and R wrapper agree exactly", {
  set.seed(42)
  # same mock data as other tests
  X0.glist <- list(matrix(rnorm(6),3,2), matrix(rnorm(6),3,2))
  X1.glist <- list(matrix(rnorm(6),3,2), matrix(rnorm(6),3,2))
  Theta.list <- list(diag(2), diag(2))
  a.i <- c(1,1); b.i <- c(0,0)
  ksi <- 0.1; gamma <- 0.1
  
  # call the R wrapper
  outR <- CordBat:::update.CorrectCoef(X0.glist, X1.glist, Theta.list,
                                       a.i, b.i, ksi, gamma, print.detail=FALSE)
  # call the C++ function directly
  outCpp <- updateCorrectCoefCpp(X0.glist, X1.glist, Theta.list,
                                 a.i, b.i, ksi, gamma)
  
  expect_equal(outR$coef.a, outCpp$coef.a)
  expect_equal(outR$coef.b, outCpp$coef.b)
})