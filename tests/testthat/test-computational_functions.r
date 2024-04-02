library(usethis)
library(testthat)

#### slice_matrix ####

test_that("Dimension of slice_matrix is correct", {
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

test_that("slice matrix works", {
  expect_equal(matrix(1:27, 9, 3) |>
                 slice_matrix(3, 3),
               array(c(1:3, 10:12, 19:21,
                        4:6, 13:15, 22:24,
                        7:9, 16:18, 25:27), c(3, 3, 3)))
})

#### invert3x3 ####
test_that("inver3x3 inverts correctly", {
  M = matrix(c(2,6,34,
               6,3,-23,
               7,-2, 0), 3, 3)
  expect_equal(M |> invert3x3() |> round(15), solve(M) |> round(15))
  M = matrix(rnorm(9), 3)
  expect_equal(M |> invert3x3() |> round(15), solve(M) |> round(15))
})

#### observation_transformation #####

