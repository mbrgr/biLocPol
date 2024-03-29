######### Libraries ###########
library(tidyverse)
library(locpol) # only used for epaK2d

# Parallel R-Code with Windows
library(parallel)
library(future.apply)


#' @title Slices a Matrix into an Array
#' @description Slices a Matrix into an Array according to the dimension of the matrix. Needed to handle the return of applying the
#'  weights.point function
#' @param M Matrix that shall be sliced
#' @param d1 rows of matrices
#' @param d2 how many matrices shall be created
#'
#' @return An Array of dimension d1 x d2 x d3, where d3 is the amount of rows of M
#' @export
#' @examples slice_matrix(M = matrix(1:27, 9, 3), d1 = 3, d2 = 3)
slice_matrix = function(M, d1, d2){
  d3 = dim(M)[2]
  if(d1 * d2 != dim(M)[1])
    stop("Slicing is not of correct dimension")
  tmp_arr = array(t(M), dim = c(d3, d1, d2))
  aperm(tmp_arr, perm = c(2, 1, 3), resize = T)
}

# epak2d = function(x, h = 1){
#
# }


