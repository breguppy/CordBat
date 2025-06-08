library(testthat)
library(CordBat)

# 1) Basic structural tests
test_that("BEgLasso returns correct structure and types", {
  set.seed(123)
  p <- 5
  n0 <- 20
  n1 <- 15
  G <- 2
  
  X0.glist <- replicate(G,
                        matrix(rnorm(n0 * p), nrow = n0, ncol = p),
                        simplify = FALSE)
  X1.glist <- lapply(X0.glist, function(mat) {
    # introduce slight batch noise
    mat + matrix(rnorm(nrow(mat) * ncol(mat), sd = 0.1), 
                 nrow = nrow(mat), ncol = ncol(mat))
  })
  
  res <- BEgLasso(
    X0.glist     = X0.glist,
    X1.glist     = X1.glist,
    penal.rho    = 0.1,
    penal.ksi    = 0.01,
    penal.gamma  = 0.01,
    eps          = 1e-3,
    print.detail = FALSE
  )
  
  # Top-level list
  expect_type(res, "list")
  expect_named(res, c("Theta", "X1.cor", "coef.a", "coef.b"))
  
  # Theta: list of p×p symmetric positive-definite matrices
  expect_type(res$Theta, "list")
  expect_length(res$Theta, G)
  for (g in seq_len(G)) {
    Theta_g <- res$Theta[[g]]
    expect_true(is.matrix(Theta_g))
    expect_equal(dim(Theta_g), c(p, p))
    # symmetry
    expect_true(max(abs(Theta_g - t(Theta_g))) < 1e-6)
    # positive definiteness
    eigs <- eigen(Theta_g, symmetric = TRUE, only.values = TRUE)$values
    expect_true(all(eigs > 0))
  }
  
  # X1.cor: list of corrected matrices with same dims as X1.glist
  expect_type(res$X1.cor, "list")
  expect_length(res$X1.cor, G)
  for (g in seq_len(G)) {
    X1c <- res$X1.cor[[g]]
    expect_equal(dim(X1c),
                 dim(X1.glist[[g]]),        # compare to the *input* dims
                 info = paste("group", g))
  }
  
  # Coefficients: numeric vectors
  expect_type(res$coef.a, "double")
  expect_length(res$coef.a, p)
  expect_true(all(res$coef.a >= 0))
  
  expect_type(res$coef.b, "double")
  expect_length(res$coef.b, p)
})

# 2) When X1 == X0, coefficients should be near (1, 0)
test_that("BEgLasso yields identity-like coefficients when no batch effect", {
  set.seed(456)
  p <- 4
  n <- 30
  X0 <- matrix(rnorm(n * p), nrow = n, ncol = p)
  # exact same data for group 1
  X1 <- X0
  
  res2 <- BEgLasso(
    X0.glist     = list(X0),
    X1.glist     = list(X0),
    penal.rho    = 0,
    penal.ksi    = 1,
    penal.gamma  = 1,
    eps          = 1e-6,
    print.detail = FALSE
  )
  expect_equal(res2$coef.a, rep(1, ncol(X0)), tolerance = 1e-6)
  expect_equal(res2$coef.b, rep(0, ncol(X0)), tolerance = 1e-6)
})

# test script for BEgLasso.R - testcases are NOT comprehensive!

# ------------------------------------------------------------------
# Helper for quickly generating toy data of matching dimensions ----
# ------------------------------------------------------------------
simulate_batches <- function(G = 2, n = 15, p = 6) {
  X0 <- lapply(seq_len(G), \(.) matrix(rnorm(n * p), n, p))
  X1 <- lapply(seq_len(G), \(.) matrix(rnorm(n * p), n, p))
  list(X0 = X0, X1 = X1, p = p, G = G)
}

test_that("print.detail toggles console output / messages", {
  set.seed(1)
  batches <- simulate_batches(G = 1, n = 10, p = 4)
  
  # Silent mode
  expect_silent(
    BEgLasso(batches$X0, batches$X1,
             penal.rho = 0.05, penal.ksi = 0.05, penal.gamma = 0.05,
             eps = 1e-2, print.detail = FALSE)
  )
})

test_that("Basic input validation — wrong shapes / types generate errors", {
  
  # Non‑list inputs ------------------------------------------------
  expect_error(
    BEgLasso(matrix(rnorm(20), 4, 5),
             list(matrix(rnorm(20), 4, 5)),
             penal.rho = 0.1, penal.ksi = 0.1, penal.gamma = 0.1,
             eps = 1e-2, print.detail = FALSE),
    regexp = "array"
  )
  
  # List length mismatch -------------------------------------------
  dat <- simulate_batches(G = 2, n = 10, p = 4)
  expect_error(
    BEgLasso(dat$X0, dat$X1[1],              # lengths differ
             penal.rho = 0.1, penal.ksi = 0.1, penal.gamma = 0.1,
             eps = 1e-2, print.detail = FALSE),
    regexp = "non-conformable"
  )
})


test_that("Output precision matrices are positive‑definite (eigenvalues > 0)", {
  set.seed(888)
  dat <- simulate_batches(G = 1, n = 20, p = 5)
  res <- BEgLasso(dat$X0, dat$X1,
                  penal.rho = 0.2, penal.ksi = 0.1, penal.gamma = 0.1,
                  eps = 1e-2, print.detail = FALSE)
  
  eig <- eigen(res$Theta[[1]], symmetric = TRUE, only.values = TRUE)$values
  expect_true(all(eig > 0))
})