library(testthat)
# skip if JGL is not installed
skip_if_not_installed("JGL")

library(CordBat)

# 1) Basic structure and types
test_that("findBestPara_JGL returns proper list structure and types", {
  set.seed(42)
  # small example: 2 groups, p=4 features, n0=10, n1=8
  p <- 4
  n0 <- 10
  n1 <- 8
  G <- 2
  
  X0.glist <- replicate(G,
                        matrix(rnorm(n0 * p), nrow = n0, ncol = p),
                        simplify = FALSE)
  X1.glist <- replicate(G,
                        matrix(rnorm(n1 * p), nrow = n1, ncol = p),
                        simplify = FALSE)
  
  lambda1.seq <- c(0.1, 1)
  lambda2.seq <- c(0.1, 1)
  
  res <- CordBat:::findBestPara_JGL(
    X0.glist    = X0.glist,
    X1.glist    = X1.glist,
    lambda1.seq = lambda1.seq,
    lambda2.seq = lambda2.seq,
    r           = 0.5
  )
  
  # top-level list
  expect_type(res, "list")
  expect_named(res, c("Theta", "best", "EBIC"))
  
  # Theta: list of G precision matrices
  expect_type(res$Theta, "list")
  expect_length(res$Theta, G)
  for (g in seq_len(G)) {
    Th <- res$Theta[[g]]
    expect_true(is.matrix(Th))
    expect_equal(dim(Th), c(p, p))
    # ensure symmetry
    expect_true(max(abs(Th - t(Th))) < 1e-6)
    # ensure positive definiteness
    ev <- eigen(Th, symmetric = TRUE, only.values = TRUE)$values
    expect_true(all(ev > 0))
  }
  
  # best: named numeric of length 2
  expect_type(res$best, "double")
  expect_length(res$best, 2)
  expect_named(res$best, c("lambda1", "lambda2"))
  
  # EBIC surface: matrix with correct dims
  expect_true(is.matrix(res$EBIC))
  expect_equal(dim(res$EBIC), c(length(lambda1.seq), length(lambda2.seq)))
})

# 2) Identity case: when X1 == X0 across two groups, no penalties should be selected
test_that("findBestPara_JGL picks zero penalties when no batch effect", {
  set.seed(123)
  p <- 3
  n <- 12
  # now use two identical groups to avoid JGL single-group bug
  X <- matrix(rnorm(n * p), nrow = n, ncol = p)
  X0.glist <- list(X, X)
  X1.glist <- list(X, X)
  
  lambda1.seq <- c(0, 1)
  lambda2.seq <- c(0, 1)
  
  res2 <- CordBat:::findBestPara_JGL(
    X0.glist    = X0.glist,
    X1.glist    = X1.glist,
    lambda1.seq = lambda1.seq,
    lambda2.seq = lambda2.seq,
    r           = 0
  )
  
  # With no structure to learn, EBIC is minimized at the most aggressive penalty
  expect_equal(res2$best["lambda1"], min(lambda1.seq))
  expect_equal(res2$best["lambda2"], min(lambda2.seq))
})
