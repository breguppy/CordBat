library(testthat)
library(CordBat)

test_that("returns 1 when there are no valid folds", {
  # N too small (no N %% k == 0 for k=2:9 with N%/%k >=10)
  expect_equal(CordBat:::selfoldforCV(19), 1)
  # divisible by 5 but quotient < 10 → still no valid folds
  expect_equal(CordBat:::selfoldforCV(25), 1)
})

test_that("picks the unique fold when exactly one divisor works", {
  # 20 %/% 2 = 10, 20 %% 2 == 0 then the only fold = 2
  expect_equal(CordBat:::selfoldforCV(20), 2)
  # 500 divisible by 5 (and by 2,4 too, but 5 is closest to 5)
  expect_equal(CordBat:::selfoldforCV(500), 5)
})

test_that("breaks ties by choosing the smaller fold when N ≤ 300", {
  # 231 is divisible only by 3 and 7, both 2 away from 5 → choose min(3,7)=3
  expect_equal(CordBat:::selfoldforCV(231), 3)
})

test_that("breaks ties by choosing the larger fold when N > 300", {
  # 399 divisible only by 3 and 7, both 2 away from 5 → choose max(3,7)=7
  expect_equal(CordBat:::selfoldforCV(399), 7)
})