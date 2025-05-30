library(testthat)
library(CordBat)

test_that("findBestPara returns correct structure and default choice on identical data", {
  skip_on_cran()
  
  set.seed(42)
  # small example: G=2 groups, p=5 features, n0=10, n1=8 samples
  G <- 2; p <- 5; n0 <- 10; n1 <- 8
  X0.glist <- lapply(seq_len(G), function(i) matrix(rnorm(n0 * p), n0, p))
  # use identical X1 so there is no real batch effect
  X1.glist <- X0.glist
  
  # run grid to pick (ksi,gamma)
  res <- CordBat:::findBestPara(
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
  skip_on_cran()
  
  set.seed(123)
  # create a clear shift in X1 so that small penal.ksi/gamma do better
  G <- 1; p <- 4; n0 <- 20; n1 <- 15
  X0 <- matrix(rnorm(n0 * p), n0, p)
  X1 <- X0 + 5  # strong offset
  res <- CordBat:::findBestPara(
    list(X0), list(X1),
    penal.rho    = 0.1,
    eps          = 1e-4,
    print.detail = FALSE
  )
  
  # we expect a smaller ksi will perform better than ksi=1
  expect_lt(res$penal.ksi, 1)
  # gamma may likewise be <1
  expect_lt(res$penal.gamma, 1)
  expect_true(is.numeric(res$MinAvedist) && is.finite(res$MinAvedist))
})
