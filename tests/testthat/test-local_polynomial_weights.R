#### weights_point ####

test_that("correct dimension", {
  p = 20
  d_grid = observation_grid(p, comp = "less")
  expect_equal(weights_point(x = c(0.2, 0.3), d_grid, h = 0.2, K = epak_2d, m = 1, del = 0) |>
                 dim(), c(1, p*(p-1)/2))
  expect_equal(weights_point(x = c(0.2, 0.3), d_grid, h = 0.2, K = epak_2d, m = 1, del = 1) |>
                 dim(), c(p*(p-1)/2, 3))
  expect_equal(weights_point(x = c(0.2, 0.3), d_grid, h = 0.2, K = epak_2d, m = 2, del = 2) |>
                 dim(), c(p*(p-1)/2, 6))
  d_grid = observation_grid(p, comp = "full")
  expect_equal(weights_point(x = c(0.2, 0.3), d_grid, h = 0.2, K = epak_2d, m = 2, del = 0) |>
                 dim(), c(1, p^2))
  expect_equal(weights_point(x = c(0.2, 0.3), d_grid, h = 0.2, K = epak_2d, m = 2, del = 1) |>
                 dim(), c(p^2, 3))
  expect_equal(weights_point(x = c(0.2, 0.3), d_grid, h = 0.2, K = epak_2d, m = 2, del = 2) |>
                 dim(), c(p^2, 6))
  d_grid = observation_grid(p, comp = "without diagonal")
  expect_equal(weights_point(x = c(0.2, 0.3), d_grid, h = 0.2, K = epak_2d, m = 1, del = 0) |>
                 dim(), c(1, p*(p-1)))
  expect_equal(weights_point(x = c(0.2, 0.3), d_grid, h = 0.2, K = epak_2d, m = 1, del = 1) |>
                 dim(), c(p*(p-1), 3))
  expect_equal(weights_point(x = c(0.2, 0.3), d_grid, h = 0.2, K = epak_2d, m = 2, del = 2) |>
                 dim(), c(p*(p-1), 6))
})

# TODO:
test_that("warning for small bandwidth", {
  p = 8
  d_grid = observation_grid(p, comp = "less")
  expect_warning(weights_point(x = c(0.2, 0.3), d_grid, h = 0.2, K = epak_2d, m = 1, del = 0))
})
