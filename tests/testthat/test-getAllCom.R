library(testthat)
library(CordBat)

test_that("getAllCom on single-feature returns list(1)", {
  X <- matrix(rnorm(10), nrow=10, ncol=1)
  comms <- CordBat:::getAllCom(X)
  expect_equal(comms, list(1L))
})

test_that("getAllCom on identical columns returns one group", {
  X <- matrix(rep(1:10, 3), nrow=10, ncol=3)  # 3 identical cols
  comms <- CordBat:::getAllCom(X)
  expect_length(comms, 1)
  expect_equal(sort(comms[[1]]), 1:3)
})

test_that("getAllCom on independent noise covers all features in small groups", {
  set.seed(42)
  X <- matrix(rnorm(100), nrow = 10, ncol = 10)
  comms <- CordBat:::getAllCom(X)
  # coverage
  expect_equal(sort(unlist(comms)), 1:10)
  # no group too large
  maxCOMSize <- min(round(nrow(X) * 0.4), 50)
  expect_true(all(lengths(comms) <= maxCOMSize))
})

test_that("getAllCom assigns every feature exactly once", {
  set.seed(2025)
  X <- matrix(rnorm(60), nrow=10, ncol=6)
  comms <- CordBat:::getAllCom(X)
  all_feats <- sort(unlist(comms))
  expect_equal(all_feats, seq_len(ncol(X)))
  expect_length(all_feats, length(unique(all_feats)))
})

test_that("getAllCom is reproducible (identical runs)", {
  set.seed(99)
  X <- matrix(rnorm(200), nrow=20, ncol=10)
  c1 <- CordBat:::getAllCom(X)
  c2 <- CordBat:::getAllCom(X)
  expect_identical(c1, c2)
})
