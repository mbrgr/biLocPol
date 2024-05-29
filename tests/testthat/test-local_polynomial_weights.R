#### weights_point ####

testthat::test_that("correct dimension", {
  p = 20
  d_grid = observation_grid(p, comp = "less")
  expect_equal(weights_point(x = c(0.2, 0.3), d_grid, h = 0.2, K = epak_2d, m = 1, del = 0) |>
                 length(), p*(p-1)/2)
  expect_equal(weights_point(x = c(0.2, 0.3), d_grid, h = 0.2, K = epak_2d, m = 1, del = 1) |>
                 dim(), c(p*(p-1)/2, 3))
  expect_equal(weights_point(x = c(0.2, 0.3), d_grid, h = 0.2, K = epak_2d, m = 2, del = 2) |>
                 dim(), c(p*(p-1)/2, 6))

  d_grid = observation_grid(p, comp = "full")
  expect_equal(weights_point(x = c(0.2, 0.3), d_grid, h = 0.2, K = epak_2d, m = 2, del = 0) |>
                 length(), p^2)
  expect_equal(weights_point(x = c(0.2, 0.3), d_grid, h = 0.2, K = epak_2d, m = 2, del = 1) |>
                 dim(), c(p^2, 3))
  expect_equal(weights_point(x = c(0.2, 0.3), d_grid, h = 0.2, K = epak_2d, m = 2, del = 2) |>
                 dim(), c(p^2, 6))
  d_grid = observation_grid(p, comp = "without diagonal")
  expect_equal(weights_point(x = c(0.2, 0.3), d_grid, h = 0.2, K = epak_2d, m = 1, del = 0) |>
                 length(), p*(p-1))
  expect_equal(weights_point(x = c(0.2, 0.3), d_grid, h = 0.2, K = epak_2d, m = 1, del = 1) |>
                 dim(), c(p*(p-1), 3))
  expect_equal(weights_point(x = c(0.2, 0.3), d_grid, h = 0.2, K = epak_2d, m = 2, del = 2) |>
                 dim(), c(p*(p-1), 6))
})

testthat::test_that("return type", {
  p = 20
  d_grid = observation_grid(p, comp = "less")
  expect_true( is.vector(weights_point(x = c(0.2, 0.3), d_grid, h = 0.2, K = epak_2d, m = 1, del = 0)))
  expect_true( is.vector(weights_point(x = c(0.2, 0.3), d_grid, h = 0.2, K = epak_2d, m = 2, del = 0)))
  expect_true( is.matrix(weights_point(x = c(0.2, 0.3), d_grid, h = 0.2, K = epak_2d, m = 1, del = 1)))
  expect_true( is.matrix(weights_point(x = c(0.2, 0.3), d_grid, h = 0.2, K = epak_2d, m = 2, del = 2)))
  expect_false( is.data.frame(weights_point(x = c(0.2, 0.3), d_grid, h = 0.2, K = epak_2d, m = 1, del = 1)))
  expect_false( is.data.frame(weights_point(x = c(0.2, 0.3), d_grid, h = 0.2, K = epak_2d, m = 2, del = 2)))
})


testthat::test_that("warning for small bandwidth", {
  p = 5
  d_grid = observation_grid(p, comp = "less")
  expect_warning(weights_point(x = c(0.2, 0.3), d_grid, h = 0.2, K = epak_2d, m = 1, del = 0))
})



#### local_polynomial_weights ####

testthat::test_that("warning for del > m", {
  p = 10
  h = 0.5
  p.eval = 5
  expect_warning(local_polynomial_weights(p, h, p.eval, m = 0, del = 1))
  expect_warning(local_polynomial_weights(p, h, p.eval, m = 0, del = 2))
  expect_warning(local_polynomial_weights(p, h, p.eval, m = 1, del = 2))
})

p = 20
h = 0.2
p.eval = 25
w = local_polynomial_weights(p, h, p.eval, m = 0)
w$weights
w2 = local_polynomial_weights(p, h, p.eval, m = 1)
w2$weights
