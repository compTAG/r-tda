context("circleUnif")

test_that("it returns required number of points", {
  pts <- circleUnif(n = 10)
  expect_equal(nrow(pts), 10)
})

test_that("points are 2d", {
  pts <- circleUnif(n = 10)
  expect_equal(ncol(pts), 2)
})

test_that("points are on circle", {
  pts <- circleUnif(n = 4, r = 2.34)
  norms <- apply(pts, MARGIN = 1, FUN = norm, type = "2")
  expect_equal(norms, rep(2.34, 4))
})
