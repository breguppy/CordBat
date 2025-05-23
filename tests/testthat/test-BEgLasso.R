library(testthat)
library(CordBat)

# Skip tests if glassoFast is not available
skip_if_not_installed("glassoFast")

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
  
  res <- CordBat:::BEgLasso(
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
  
  res2 <- CordBat:::BEgLasso(
    X0.glist     = list(X0),
    X1.glist     = list(X0),
    penal.rho    = 0,
    penal.ksi    = 0,
    penal.gamma  = 0,
    eps          = 1e-6,
    print.detail = FALSE
  )
  expect_equal(res2$coef.a, rep(1, ncol(X0)), tolerance = 1e-6)
  expect_equal(res2$coef.b, rep(0, ncol(X0)), tolerance = 1e-6)
})