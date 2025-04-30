library(testthat)
library(CordBat)

test_that("ComtyDet leaves small communities unchanged", {
  # G = 4×4 identity => no edges, everything is size <= 2
  G <- diag(1, 4)
  InputCOM <- list(1:2, 3:4)
  out <- CordBat:::ComtyDet(G, InputCOM, minCOMSize = 2)
  expect_equal(out, InputCOM)
})

test_that("ComtyDet splits a clear two-block structure", {
  # two perfect blocks and no cross-block edges
  G <- matrix(0,4,4)
  G[1,2] <- G[2,1] <- 1
  G[3,4] <- G[4,3] <- 1
  
  InputCOM <- list(1:4)
  out <- CordBat:::ComtyDet(G, InputCOM, minCOMSize = 1)
  
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
  out <- CordBat:::ComtyDet(G, InputCOM, minCOMSize = 2)
  
  all_in <- sort(unlist(InputCOM))
  all_out <- sort(unlist(out))
  expect_equal(all_out, all_in)
  expect_length(all_out, length(unique(all_out)))
})