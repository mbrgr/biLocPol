#### Computational functions ####
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
slice_matrix = function(M, d1, d2) {
  d3 = dim(M)[2]
  if(d1 * d2 != dim(M)[1])
    stop("Slicing is not of correct dimension")
  tmp_arr = array(t(M), dim = c(d3, d1, d2))
  aperm(tmp_arr, perm = c(2, 1, 3), resize = T)
}


#' Inverts a 3x3 matrix
#' @description Inverts a 3x3 matrix more quickly than solve.
#'
#' @param M Matrix of Dimension 3x3
#' @param eps If the the determinant of M is smaller than this value the function returns an error (unstable matrix inversion).
#'
#' @return The inverse 3x3 matrix.
#' @export
#'
#' @examples
#' M = matrix(c(2,3,5,2,5,3,5,3,2), 3, 3)
#' M_inv = invert3x3(M)
#' round(M %*% M_inv, 15)
invert3x3 = function(M, eps = 10^(-6)) {
  stopifnot("M is not a matrix" = is.matrix(M))
  stopifnot("M is not a 3x3 matrix" = all(dim(M) == c(3, 3)) )

  B = matrix(0, 3, 3)

  N_11 = (M[2, 3] * M[3, 2] - M[2, 2] * M[3, 3])
  N_21 = (M[2, 1] * M[3, 3] - M[2, 3] * M[3, 1])
  N_31 = (M[2, 1] * M[3, 2] - M[2, 2] * M[3, 1])
  N_12 = (M[1, 3] * M[3, 2] - M[1, 2] * M[3, 3])
  N_22 = (M[1, 3] * M[3, 1] - M[1, 1] * M[3, 3])
  N_32 = (M[1, 1] * M[3, 2] - M[1, 2] * M[3, 1])
  N_13 = (M[1, 3] * M[2, 2] - M[1, 2] * M[2, 3])
  N_23 = (M[1, 1] * M[2, 3] - M[1, 3] * M[2, 1])
  N_33 = (M[1, 1] * M[2, 2] - M[1, 2] * M[2, 1])
  det_M = M[1, 1] * N_11 - M[2, 1] * N_12 + M[3, 1] * N_13

  if(abs(det_M) < eps){
    print("unstable matrix inversion")
  }
  minus_det_M = -det_M
  det_M_inv = 1/det_M
  minus_det_M_inv = 1/minus_det_M

  B[1, 1] = N_11 * det_M_inv
  B[2, 1] = N_21 * det_M_inv
  B[3, 1] = N_31 * minus_det_M_inv
  B[1, 2] = N_12 * minus_det_M_inv
  B[2, 2] = N_22 * det_M_inv
  B[3, 2] = N_32 * det_M_inv
  B[1, 3] = N_13 * det_M_inv
  B[2, 3] = N_23 * det_M_inv
  B[3, 3] = N_33 * minus_det_M_inv
  return(B)
}


#### Kernel and U ####

#' Two-Dimensional Epanechnikov Kernel
#' @description Two dimensional Epanechnikov Kernel with optional bandwidth parameter.
#'
#' @param x Two dimensional vector at which the kernel shall be evaluated.
#' @param h Bandwidth parameter for nonparametric regression. Default is h = 1.
#'
#' @return Real value between 0 and 1.
#' @export
#'
#' @examples
#' epak_2d(c(0.2, 0.13), 0.4)
epak_2d = function(x, h = 1) {
  x = x/h
  (2/pi) * (1 - x[1]^2 - x[2]^2) * ((x[1]^2 + x[2]^2) < 1)
}


#' Vector U for bivariate local polynomial estimator
#'
#' @param x two dimensional vector of evaluation.
#' @param h bandwidth parameter.
#' @param m `m = 0`refers to Nadaraya-Watson, `m = 1` to the local linear estimator and `m = 2` to local quadratic
#'
#' @return one, three or six dimensional vector (depending on `m`)
#' @export
#'
#' @examples
#' U(rep(0.5, 6), h = 0.25, m = 2)
U = function(x, h, m = 1){
  x = x/h
  if (m == 0) {
    return(1)
  } else if (m == 1) {
    return(c(1, x[1], x[2]))
  } else if (m == 2) {
    return(c(1, x[1], x[2], x[1]^2/2, x[1]*x[2], x[2]^2/2))
  } else {
    stop("m not implemented")
  }
}



#### Creating and Transforming Oberservations and Grids ####

#' Transforms observations according to Berger/Holzmann (2024)
#'
#' @param Y Observations in Form of a n x p matrix
#' @param grid.type which kind of observation transformation (`less`, `without diagonal`, `full`)
#'
#' @return Vector with transformed obseravtions of dimension `p(p-1)/2` or `p(p-1)` or `p^2`
#' @export
#'
#' @examples
#' 0
observation_transformation = function(Y, grid.type = "less"){
  p = length(Y[1, ])
  n = length(Y[, 1])
  Y.means = colMeans(Y)                         # pointwise mean
  Y.2 = apply(Y, 1, tcrossprod)                 # Dim (p^2 x n) (order of the col's are as.vector of the matrices)
  Z = Y.2 - as.vector(tcrossprod(Y.means))      # vector gets subtracted of all columns. Dim (p^2 x n)
  Z.mean = rowSums(Z)/(n-1)                     # Dim p^2. Reduction that is only possible in the case of synchronous observations.
  if (grid.type == "less") {
    M2 = as.vector(upper.tri(matrix(0, p, p)))    # corresponds to entries of: observation.grid(p, "less)
  } else if (grid.type == "without diagonal") {
    M2 = as.vector(!diag(T, p, p))
  } else if (grid.type == "full") {
    M2 = T
  } else {stop("grid type not implemented")}
  return(Z.mean[M2])                                    # in
}

#' Creates two dimensional observation grids
#' @description Creates two dimensional observation grids.
#'
#' @param p Amount of design points.
#' @param x Needed if p is not given. p-dimensional vector with design points on a real interval. If x
#' is not given x is created with p equidistant design points.
#' @param comp
#' Indicates the structure of the grid:
#'
#' -- `less`: The first coordinate from the grid is always strictly smaller then the second, `gtr`: vice versa.
#'
#' -- `lesseq`:  The first coordinate from the grid is always smaller or equal then the second, `gtreq`: vice versa.
#'
#' -- `without diagonal`: Elements with equal values are excluded from the grid.
#'
#' -- `full`: No element is excluded.
#'
#' @return Matrix of dimension q x 2 where q depends on comp. Note the rownames of the matrix refer to the initial
#' numbering before excluding certain values as indicated by comp.
#' @export
#'
#' @examples
#' observation_grid(10)
observation_grid = function(p = NULL, x = NULL, comp = "less") {
  if(is.null(x) & is.null(p))
    stop("p or grid needs to be supplied")
  else if(is.null(x))
    x = (1:p - 1/2)/p
  x.grid = expand.grid(x, x)

  if (comp == "less")
    b = x.grid[,1] < x.grid[,2]
  if (comp == "lesseq")
    b = x.grid[,1] <= x.grid[,2]
  if (comp == "gtr")
    b = x.grid[,1] > x.grid[,2]
  if (comp == "gtreq")
    b = x.grid[,1] >= x.grid[,2]
  if (comp == "without diagonal")
    b = !(x.grid[,1] == x.grid[,2])
  if (comp == "full")
    b = T
  x.grid[b, ]
}

