library(testthat)
library(CordBat)   # load your package

test_that("soft returns zero when abs(x) <= lambda", {
  # scalars
  expect_equal(CordBat:::soft(0,   1), 0)
  expect_equal(CordBat:::soft(0.5, 1), 0)
  expect_equal(CordBat:::soft(-0.5,1), 0)
  
  # vector
  expect_equal(CordBat:::soft(c(-0.2, 0.2), 0.5),
               c(0, 0))
})

test_that("soft shrinks magnitude by lambda when abs(x) > labmda", {
  expect_equal(CordBat:::soft(3, 1),   2)
  expect_equal(CordBat:::soft(-3, 1), -2)
})

test_that("soft is vectorized", {
  x <- c(-3, -1,  0,  1,  3)
  lambda <- 1
  expect_equal(CordBat:::soft(x, lambda),
               c(-2,  0,  0,  0,  2))
})

test_that("lambda = 0 leaves x unchanged", {
  expect_equal(CordBat:::soft(5,       0),  5)
  expect_equal(CordBat:::soft(-2:2,   0), -2:2)
})

test_that("lambda > max(abs(x)) maps everything to zero", {
  x <- c(-1, 0, 1)
  expect_equal(CordBat:::soft(x, 2),
               c(0, 0, 0))
})