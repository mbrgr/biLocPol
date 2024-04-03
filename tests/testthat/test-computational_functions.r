#### slice_matrix ####

testthat::test_that("Dimension of slice_matrix is correct", {
 expect_equal(matrix(1:36, 9, 4) |>
                slice_matrix(3, 3) |> dim(),
              c(3, 4, 3))
  expect_equal(matrix(1:36, 12, 3) |>
                 slice_matrix(3, 4) |> dim(),
               c(3, 3, 4))
  expect_equal(matrix(1:36, 12, 3) |>
                 slice_matrix(4, 3) |> dim(),
               c(4, 3, 3))
})

testthat::test_that("slice matrix works", {
  expect_equal(matrix(1:27, 9, 3) |>
                 slice_matrix(3, 3),
               array(c(1:3, 10:12, 19:21,
                        4:6, 13:15, 22:24,
                        7:9, 16:18, 25:27), c(3, 3, 3)))
})

#### invert3x3 ####
testthat::test_that("inver3x3 inverts correctly", {
  M = matrix(c(2,  6,  34,
               6,  3, -23,
               7, -2,   0), 3, 3)
  expect_equal(invert3x3(M), solve(M))
  M = matrix(rnorm(9), 3)
  expect_equal(invert3x3(M), solve(M))
})

testthat::test_that("expect warning", {
  expect_warning(matrix(1, 3, 3) |> invert3x3())
})

#### observation_transformation #####
testthat::test_that("correct dimension", {
  n = 20
  p = 40
  Y = FDA_observation(n, (1:p - 0.5)/p)
  expect_length(Y |> observation_transformation(grid.type = "less"), p*(p-1)/2)
  expect_length(Y |> observation_transformation(grid.type = "without diagonal"), p*(p-1))
  expect_length(Y |> observation_transformation(grid.type = "full"), p^2)
})

testthat::test_that("correct values", {
  M = matrix(1:12, 3, 4)
  expect_equal(M |> observation_transformation(), rep(1, 6))
})

#### observation_grid ####
testthat::test_that("correct dimension", {
  p = 20
  expect_equal( observation_grid(p, comp = "less") |> dim(), c(p*(p-1)/2, 2))
  expect_equal( observation_grid(p, comp = "lesseq") |> dim(), c(p*(p+1)/2, 2))
  expect_equal( observation_grid(p, comp = "without diagonal") |> dim(), c(p*(p-1), 2))
  expect_equal( observation_grid(p, comp = "gtr") |> dim(), c(p*(p-1)/2, 2))
  expect_equal( observation_grid(p, comp = "gtreq") |> dim(), c(p*(p+1)/2, 2))
  expect_equal( observation_grid(p, comp = "full") |> dim(), c(p^2, 2))
})
