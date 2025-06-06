library(testthat)
library(CordBat)

test_that("getAllCom on single-feature returns list(1)", {
  X <- matrix(rnorm(10), nrow=10, ncol=1)
  comms <- getAllCom(X)
  expect_equal(comms, list(1L))
})

test_that("getAllCom on identical columns returns one group", {
  X <- matrix(rep(1:10, 3), nrow=10, ncol=3)  # 3 identical cols
  comms <- getAllCom(X)
  expect_length(comms, 1)
  expect_equal(sort(comms[[1]]), 1:3)
})

test_that("getAllCom on independent noise covers all features in small groups", {
  set.seed(42)
  X <- matrix(rnorm(100), nrow = 10, ncol = 10)
  comms <- getAllCom(X)
  # coverage
  expect_equal(sort(unlist(comms)), 1:10)
  # no group too large
  maxCOMSize <- min(round(nrow(X) * 0.4), 50)
  expect_true(all(lengths(comms) <= maxCOMSize))
})

test_that("getAllCom assigns every feature exactly once", {
  set.seed(2025)
  X <- matrix(rnorm(60), nrow=10, ncol=6)
  comms <- getAllCom(X)
  all_feats <- sort(unlist(comms))
  expect_equal(all_feats, seq_len(ncol(X)))
  expect_length(all_feats, length(unique(all_feats)))
})

test_that("getAllCom is reproducible (identical runs)", {
  set.seed(99)
  X <- matrix(rnorm(200), nrow=20, ncol=10)
  c1 <- getAllCom(X)
  c2 <- getAllCom(X)
  expect_identical(c1, c2)
})

test_that("ComtyDet leaves small communities unchanged", {
  # G = 4×4 identity => no edges, everything is size <= 2
  G <- diag(1, 4)
  InputCOM <- list(1:2, 3:4)
  out <- ComtyDet(G, InputCOM, minCOMSize = 2)
  expect_equal(out, InputCOM)
})

test_that("ComtyDet splits a clear two-block structure", {
  # two perfect blocks and no cross-block edges
  G <- matrix(0,4,4)
  G[1,2] <- G[2,1] <- 1
  G[3,4] <- G[4,3] <- 1
  
  InputCOM <- list(1:4)
  out <- ComtyDet(G, InputCOM, minCOMSize = 1)
  
  # we should get exactly 2 communities
  expect_length(out, 2)
  
  # and one must be {1,2}, the other {3,4} (in any order)
  comm_sets <- lapply(out, sort)
  expect_true(any(identical(comm_sets[[1]], c(1L,2L)),
                  identical(comm_sets[[2]], c(1L,2L))))
  expect_true(any(identical(comm_sets[[1]], c(3L,4L)),
                  identical(comm_sets[[2]], c(3L,4L))))
})

test_that("ComtyDet output covers all input indices exactly once", {
  G <- diag(1, 5)   # trivial graph
  InputCOM <- list(1:3, 4:5)
  out <- ComtyDet(G, InputCOM, minCOMSize = 2)
  
  all_in <- sort(unlist(InputCOM))
  all_out <- sort(unlist(out))
  expect_equal(all_out, all_in)
  expect_length(all_out, length(unique(all_out)))
})

test_that("ComtyDet returns identical output for a fixed graph", {
  set.seed(42)  
  
  # Create a toy symmetric matrix G (e.g., correlation or similarity matrix)
  n <- 10
  G <- matrix(runif(n^2), n, n)
  G <- (G + t(G)) / 2  # make symmetric
  diag(G) <- 1         # set diagonal to 1 for simplicity
  
  InputCOM <- list(1:5, 6:10)
  minCOMSize <- 2
  
  res1 <- ComtyDet(G, InputCOM, minCOMSize)
  res2 <- ComtyDet(G, InputCOM, minCOMSize)
  
  # Compare results: same length and same content
  expect_equal(length(res1), length(res2))
  expect_equal(res1, res2)
})