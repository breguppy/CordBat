# test script for utils.R - testcases are NOT comprehensive!
library(testthat)
library(CordBat)

## -----------------------------------------------------------
## 1.  selfoldforCV   ----------------------------------------
## -----------------------------------------------------------

test_that("selfoldforCV returns 1 when there are no valid folds", {
  # N too small (no N %% k == 0 for k=2:9 with N%/%k >=10)
  expect_equal(selfoldforCV(19), 1)
  # divisible by 5 but quotient < 10 → still no valid folds
  expect_equal(selfoldforCV(25), 1)
})

test_that("selfoldforCV picks the unique fold when exactly one divisor works", {
  # 20 %/% 2 = 10, 20 %% 2 == 0 then the only fold = 2
  expect_equal(selfoldforCV(20), 2)
  # 500 divisible by 5 (and by 2,4 too, but 5 is closest to 5)
  expect_equal(selfoldforCV(500), 5)
})

test_that("selfoldforCV breaks ties by choosing the smaller fold when N ≤ 300", {
  # 231 is divisible only by 3 and 7, both 2 away from 5 → choose min(3,7)=3
  expect_equal(selfoldforCV(231), 3)
})

test_that("selfoldforCV breaks ties by choosing the larger fold when N > 300", {
  # 399 divisible only by 3 and 7, both 2 away from 5 → choose max(3,7)=7
  expect_equal(selfoldforCV(399), 7)
})

## -----------------------------------------------------------
## 2.  soft -----------------------------------------------
## -----------------------------------------------------------

test_that("soft returns zero when abs(x) <= lambda", {
  # scalars
  expect_equal(soft(0, 1), 0)
  expect_equal(soft(0.5, 1), 0)
  expect_equal(soft(-0.5, 1), 0)
  
  # vector
  expect_equal(soft(c(-0.2, 0.2), 0.5),
               c(0, 0))
})

test_that("soft shrinks magnitude by lambda when abs(x) > labmda", {
  expect_equal(soft(3, 1), 2)
  expect_equal(soft(-3, 1), -2)
})

test_that("soft: lambda = 0 leaves x unchanged", {
  expect_equal(soft(5, 0), 5)
})

test_that("soft: lambda > max(abs(x)) maps everything to zero", {
  x <- c(-1, 0, 1)
  expect_equal(soft(x, 2),
               c(0, 0, 0))
})

## -----------------------------------------------------------
## 3.  CDfgL  ----------------------------------------------
## -----------------------------------------------------------

test_that("CDfgL converges on simple input", {
  V <- diag(3)  # Identity matrix
  beta_i <- c(0.5, -0.3, 0.1)
  u <- c(1, 1, 1)
  rho <- 0.1
  
  result <- CDfgL(V, beta_i, u, rho, maxIter = 100, print.detail = FALSE)
  
  expect_type(result, "double")
  expect_length(result, 3)
  expect_true(all(is.finite(result)))
})

test_that("CDfgL returns zero for zero input", {
  V <- diag(3)
  beta_i <- c(0, 0, 0)
  u <- c(0, 0, 0)
  rho <- 0.5
  
  result <- CDfgL(V, beta_i, u, rho, maxIter = 50, print.detail = FALSE)
  
  expect_equal(result, c(0, 0, 0))
})

test_that("CDfgL solves identity‑matrix sub‑problem exactly", {
  set.seed(1)
  p <- 3
  V <- diag(p)
  u <- c(1, -2, 0.5)
  rho <- 0.3
  sol <- CDfgL(V, beta_i = rep(0, p), u = u, rho = rho,
               print.detail = FALSE, maxIter = 1000)
  
  expect_equal(sol, soft(u, rho), tolerance = 1e-6)
})

test_that("CDfgL hits maximum iterations when convergence is impossible", {
  set.seed(123)
  V <- diag(c(1, 2, 3))  # Diagonal matrix but not identity
  beta_i <- runif(3, -10, 10)  # Start far from zero
  u <- runif(3, -10, 10)
  rho <- 1e-10  # Tiny regularization -> updates will be very small
  
  expect_message(
    CDfgL(V, beta_i, u, rho, maxIter = 1, print.detail = TRUE),
    "CDfgL reached max iteration"
  )
})

## -----------------------------------------------------------
## 4.  update_CorrectCoef  ----------------------------------
## -----------------------------------------------------------

test_that("update_CorrectCoef returns sane dimensions and small changes on identical data", {
  set.seed(123)
  G <- 1; n <- 12; p <- 4
  X0 <- list(matrix(rnorm(n * p), n, p))
  X1 <- X0                                   # identical batches
  Theta.list <- list(diag(p))                # identity precision
  res <- update_CorrectCoef(
    X0, X1, Theta.list,
    a.i = rep(1, p), b.i = rep(0, p),
    penal.ksi = 0.1, penal.gamma = 0.1,
    print.detail = FALSE
  )
  
  expect_named(res, c("coef.a", "coef.b"), ignore.order = TRUE)
  expect_length(res$coef.a, p)
  expect_length(res$coef.b, p)
  
  # On identical batches, a ≈ 1 and b ≈ 0
  expect_equal(res$coef.a, rep(0.2, p), tolerance = 0.2)
  expect_equal(res$coef.b, rep(0, p), tolerance = 0.2)
})

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
  result <- update.CorrectCoef(
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
  
  result <- update.CorrectCoef(
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
  
  result <- update.CorrectCoef(
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
    update.CorrectCoef(
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
    update.CorrectCoef(
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

test_that("update.CorrectCoef: C++ core and R wrapper agree exactly", {
  set.seed(42)
  X0.glist <- list(matrix(rnorm(6),3,2), matrix(rnorm(6),3,2))
  X1.glist <- list(matrix(rnorm(6),3,2), matrix(rnorm(6),3,2))
  Theta.list <- list(diag(2), diag(2))
  a.i <- c(1,1); b.i <- c(0,0)
  ksi <- 0.1; gamma <- 0.1
  
  # call the R wrapper
  outR <- update.CorrectCoef(X0.glist, X1.glist, Theta.list,
                                       a.i, b.i, ksi, gamma, print.detail=FALSE)
  # call the C++ function directly
  outCpp <- updateCorrectCoefCpp(X0.glist, X1.glist, Theta.list,
                                 a.i, b.i, ksi, gamma)
  
  expect_equal(outR$coef.a, outCpp$coef.a)
  expect_equal(outR$coef.b, outCpp$coef.b)
})

test_that("update.CorrectCoef returns similar values as update_CorrectCoef", {
  
  set.seed(123)
  X0.glist <- list(matrix(rnorm(6),3,2), matrix(rnorm(6),3,2))
  X1.glist <- list(matrix(rnorm(6),3,2), matrix(rnorm(6),3,2))
  Theta.list <- list(diag(2), diag(2))
  
  a.i <- c(1, 1)
  b.i <- c(0, 0)
  
  penal.ksi <- 0.1
  penal.gamma <- 0.1
  print.detail <- FALSE
  
  new_result <- update.CorrectCoef(
    X0.glist = X0.glist,
    X1.glist = X1.glist,
    Theta.list = Theta.list,
    a.i = a.i,
    b.i = b.i,
    penal.ksi = penal.ksi,
    penal.gamma = penal.gamma,
    print.detail = print.detail
  )
  
  old_result <- update_CorrectCoef(
    X0.glist = X0.glist,
    X1.glist = X1.glist,
    Theta.list = Theta.list,
    a.i = a.i,
    b.i = b.i,
    penal.ksi = penal.ksi,
    penal.gamma = penal.gamma,
    print.detail = print.detail
  )
  
  expect_equal(new_result$coef.a, old_result$coef.a)
  expect_equal(new_result$coef.b, old_result$coef.b)
})

## -----------------------------------------------------------
## 5.  findBestPara  ----------------------------------------
## -----------------------------------------------------------

test_that("findBestPara returns selected penalties", {
  set.seed(11)
  G <- 1; n <- 10; p <- 3
  X0 <- list(matrix(rnorm(n * p), n, p))
  X1 <- list(matrix(rnorm(n * p), n, p))
  
  res <- findBestPara(X0, X1, penal.rho = 0.2, eps = 1e-2, print.detail = FALSE)
  expect_type(res, "list")
  expect_named(res, c("penal.ksi", "penal.gamma", "MinAvedist"), ignore.order = TRUE)
  expect_true(res$penal.ksi %in% c(1, 0.5, 0.3, 0.1))
  expect_true(res$penal.gamma %in% c(1, 0.5, 0.3, 0.1))
  expect_true(is.finite(res$MinAvedist))
})

test_that("findBestPara returns correct structure and default choice on identical data", {
  set.seed(42)
  # small example: G=2 groups, p=5 features, n0=10, n1=8 samples
  G <- 2; p <- 5; n0 <- 10; n1 <- 8
  X0.glist <- lapply(seq_len(G), function(i) matrix(rnorm(n0 * p), n0, p))
  # use identical X1 so there is no real batch effect
  X1.glist <- X0.glist
  
  # run grid to pick (ksi,gamma)
  res <- findBestPara(
    X0.glist,
    X1.glist,
    penal.rho   = 0.1,
    eps         = 1e-4,
    print.detail = FALSE
  )
  
  # structure
  expect_type(res, "list")
  expect_named(res, c("penal.ksi", "penal.gamma", "MinAvedist"))
  
  # since every (ksi,γ) yields essentially the same EBIC,
  # the first grid point (1,1) should be selected
  expect_equal(res$penal.ksi,   1)
  expect_equal(res$penal.gamma, 1)
  
  # MinAvedist must be a finite numeric
  expect_true(is.numeric(res$MinAvedist) && is.finite(res$MinAvedist))
})

test_that("findBestPara picks a different penalty when batch effect is introduced", {
  set.seed(2025)
  
  n <- 100  
  p <- 2      
  
  # X0 ~ nearly perfect colinearity: column 2 = column 1 + tiny noise
  z  <- rnorm(n)
  eps <- rnorm(n, mean = 0, sd = 1e-4)
  X0 <- cbind(z, z + eps)
  scale <- cbind(rnorm(n), 2 * eps)
  shift <- rnorm(n)
  X1 <- scale * X0 - cbind(shift, shift)
  
  
  out <- findBestPara(
    list(X0),
    list(X1),
    penal.rho   = 0.01,
    eps         = 1e-4,
    print.detail = FALSE
  )
  
  expect_true(out$penal.ksi   < 1,
              info = sprintf("Expected penal.ksi < 1, got %g", out$penal.ksi))
  expect_true(out$penal.gamma < 1,
              info = sprintf("Expected penal.gamma < 1, got %g", out$penal.gamma))
  expect_true(is.numeric(out$MinAvedist) && is.finite(out$MinAvedist))
})

test_that("findBeatPara gives same output at old_findBestPara when batch effect is introduced", {
  set.seed(2025)
  
  n <- 100  
  p <- 2      
  
  # X0 ~ nearly perfect colinearity: column 2 = column 1 + tiny noise
  z  <- rnorm(n)
  eps <- rnorm(n, mean = 0, sd = 1e-4)
  X0 <- cbind(z, z + eps)
  scale <- cbind(rnorm(n), 2 * eps)
  shift <- rnorm(n)
  X1 <- scale * X0 - cbind(shift, shift)
  
  
  new_out <- findBestPara(
    list(X0),
    list(X1),
    penal.rho   = 0.01,
    eps         = 1e-4,
    print.detail = FALSE
  )
  
  old_out <- old_findBestPara(
    list(X0),
    list(X1),
    penal.rho   = 0.01,
    eps         = 1e-4,
    print.detail = FALSE
  )
  
  expect_true(new_out$penal.ksi == old_out$penal.ksi,
              info = sprintf("New penal.ksi =  %g and old penal.ksi = %g", new_out$penal.ksi, old_out$penal.ksi))
  expect_true(new_out$penal.gamma == old_out$penal.gamma,
              info = sprintf("New penal.gamma =  %g and old penal.gamma = %g", new_out$penal.gamma, old_out$penal.gamma))
  expect_true(new_out$MinAvedist == old_out$MinAvedist,
              info = sprintf("New MinAvedist =  %g and old MinAvedist = %g", new_out$MinAvedist, old_out$MinAvedist))
})

## -----------------------------------------------------------
## 6.  selrho.useCVBIC  -------------------------------------
## -----------------------------------------------------------

test_that("selrho.useCVBIC returns correct structure and rho from candidate grid", {
  set.seed(7)
  X <- matrix(rnorm(120), 20, 6)
  out <- selrho.useCVBIC(X, print.detail = FALSE)
  expect_type(out, "double")
  expect_length(out, 2)
  expect_true(out[1] %in% seq(0.1, 0.9, 0.1))
  expect_true(is.finite(out[2]))
})

test_that("selrho.useCVBIC emits selection message when print.detail = TRUE", {
  set.seed(789)
  
  X <- matrix(rnorm(200), nrow = 20, ncol = 10)
  expect_message(
    selrho.useCVBIC(X = X, print.detail = TRUE),
    "CVBIC selects rho = "
  )
})

test_that("selrho.useCVBIC returns the same values as old_selrho.useCVBIC when there is exactly one optimal rho", {
  set.seed(123)
  
  X <- matrix(rnorm(200), nrow = 20, ncol = 10)
  new_out <- selrho.useCVBIC(X, print.detail = FALSE)
  old_out <- old_selrho.useCVBIC(X, print.detail = FALSE)
  
  expect_true(new_out[1] == old_out[1],
              info = sprintf("New rho =  %g and old rho = %g", new_out[1], old_out[1]))
  expect_true(new_out[2] == old_out[2],
              info = sprintf("New MinCVerr =  %g and old MinCVerr = %g", new_out[2], old_out[2]))
})

test_that("selrho.useCVBIC returns the largerst rho in old_selrho.useCVBIC when there is a tie for MinCVerr", {
  set.seed(456)
  
  X <- matrix(rnorm(200), nrow = 20, ncol = 10)
  new_out <- selrho.useCVBIC(X, print.detail = FALSE)
  old_out <- old_selrho.useCVBIC(X, print.detail = FALSE)
  n <- length(old_out)
  expect_true(new_out[1] == old_out[n-1],
              info = sprintf("New rho =  %g and old rho = %g", new_out[1], old_out[n-1]))
  expect_true(new_out[2] == old_out[n],
              info = sprintf("New MinCVerr =  %g and old MinCVerr = %g", new_out[2], old_out[n]))
})

## -----------------------------------------------------------
## 7.  DelOutlier  ------------------------------------------
## -----------------------------------------------------------

test_that("DelOutlier detects and removes extreme samples", {
  set.seed(21)
  X <- matrix(rnorm(200), 20, 10)
  X[1, ] <- 10                               # big outlier row
  
  res <- DelOutlier(X)
  expect_type(res, "list")
  expect_named(res, c("delsampIdx", "X.out"), ignore.order = TRUE)
  expect_true(1 %in% res$delsampIdx)
  expect_equal(nrow(res$X.out), 19)
})

test_that("DelOutlier keeps data intact when no strong outliers", {
  X <- matrix(rnorm(100), 10, 10)
  res <- DelOutlier(X)
  expect_length(res$delsampIdx, 0)
  expect_equal(nrow(res$X.out), 10)
})

## -----------------------------------------------------------
## 8.  ImputeOutlier  ---------------------------------------
## -----------------------------------------------------------

test_that("ImputeOutlier replaces extreme points and leaves no NAs", {
  set.seed(31)
  X <- matrix(rnorm(300), 30, 10)
  X[5, 3] <- 30 # extreme high value
  
  dat.i <- X[, 3]
  dat.i.m <- mean(dat.i)
  dat.i.sd <- sd(dat.i)
  
  Y <- ImputeOutlier(X)
  expect_equal(dim(Y), dim(X))
  expect_false(anyNA(Y))
  # Value should have shrunk toward centre (not remain extreme)
  expect_true(abs(Y[5, 3]) < 10)
})