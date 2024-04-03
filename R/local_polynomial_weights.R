#### Calculation of the weights ####
#' @title Calculate the weights of bivariate local polynomial estimator at a single point
#'
#' @description function that gets called from local.polynomial.weights
#'  to calculate the weights at a specific point. Therefore the function is not
#'  meant to be used singularly.
#'
#' @param x design point as two dimensional vector at which the weights are calculated
#' @param x.design.grid grid of design points for calculation of the weights
#' @param h bandwidth parameter that needs to be chosen
#' @param K kernel function - default is two dimensional Epanechnikov kernel
#' @param m 1 refers to local linear, 2 to local quadratic
#' @param del partial derivatives of which order shall be calculated?
#'
#' @return returns a matrix with one to six columns where the first column
#'  is the estimation of the actual function, the second is the partial derivative
#'  in the first direction (x), the third in direction 2 (y), the fourth xx, the
#'  fifth xy and the sixth in direction yy. In case `del = 0` its an row vector and other wise a matrix.
#'
#' @export
#'
#' @examples
#' weights_point(c(0.2, 0.3), observation_grid(15, comp = "less"), 0.2, del = 1)
weights_point = function(x, x.design.grid, h, K = epak_2d, m = 1, del = 0){
 L = apply(x.design.grid, 1, function(z){tcrossprod(U(z - x, h, m = m), U(z-x, h, m = m)) * K((z-x)/h)})
 if (m == 2) {
   B = matrix(rowSums(L), 6, 6) # This differs from local linear estimator
   B.inv = solve(B) # TODO: can this be sped up????
   U_del = switch(del + 1,
                  c(T, rep(F, 5)),
                  rep(c(T, F), each = 3),
                  rep(T, 6))
 } else { # m = 1
   B = matrix(rowSums(L), 3, 3)
   B.inv = invert3x3(B)
   U_del = switch(del + 1, c(T,F,F), rep(T, 3))
 }

 by_h = 1/h
 u = switch(del+1,
             1,
             c(1, by_h, by_h),
             c(1, by_h, by_h, by_h^2, by_h^2, by_h^2))

  t(u * apply(x.design.grid, 1, function(z){
    crossprod(B.inv[,U_del], U(z-x, h, m = m) * K((z-x)/h))}))
}



# p amount of Design points -> generates equidistant grid
# p.eval are the amount of evaluation points -> best to be chosen high, till now
# only p.eval[,1] <= p.eval[,2] gets generated
# Returns array Dim p*(p-1)/2 x q x d (where q = amount of evaluation points, and d the amount of derivatives demanded)
# Result: rows are sort like observation.grid(p, "less")'s first column
# columns are sort like observation.grid(p.eval, "lesseq")' first column
# result can be evaluated by observations transformed with "observation.transform" by calling tcrossprod(t(W), Z)



#' Calculated the weights of the linear local polynomial estimator
#' @description according to paper Berger/Holzmann (2024)
#'
#' @param p amount of design (observation) points on one axis
#' @param h bandwidth parameter
#' @param p.eval amount of evaluation points (can be NULL if `eval.type = diagonal`)
#' @param parallel Paralization via the future package
#' @param m `m = 0`refers to Nadaraya-Watson, `m = 1` to the local linear estimator and `m = 2` to local quadratic
#' @param del Calculation of derivatives
#' @param x.design.grid Can be used for custom designs
#' @param grid.type which kind of grid (`less`, `without diagonal`, `full`)
#' @param eval.type at all evaluation points or only at the diagonal (`full`, `diagonal`)
#' @param ... further arguments
#'
#' @return TODO: write something usefull -> can we create an object here?
#' @importFrom future.apply future_apply
#' @importFrom parallel detectCores makeCluster stopCluster
#' @importFrom future plan
#' @export
#'
#' @examples
#' local_polynomial_weights(20, 0.3, 25, m = 2, del = 0)
local_polynomial_weights = function(p, h, p.eval, parallel = F, m = 1,
                                    del = 0, x.design.grid = NULL,
                                    grid.type = "less", eval.type = "full", ...){

  if( !(grid.type %in% c("less", "without diagonal")) ){
    stop("grid type is not feasible - choose -less- oder -without diagonal-.")
  }


  if(!is.null(x.design.grid) & is.vector(x.design.grid)){
    x.design.grid = observation_grid(x = x.design.grid, comp = grid.type)
  }else{
    x.design.grid = observation_grid(p = p, comp = grid.type)
  }

  if( !(eval.type %in% c("full", "diagonal")) ){
    stop("evaluation is only possible for lower diagonal -full- oder only on the -diagonal-.")
  }
  x.eval.grid = switch(eval.type,
                       full = observation_grid(p = p.eval, comp = "lesseq"),
                       diagonal = matrix( (1:p.eval - 0.5)/p.eval, p.eval, 2)   )

  if(del > m){m = 2; del = 2}
  if(parallel){
    cl = parallel::makeCluster(parallel::detectCores( ) - 1)
    future::plan(future::cluster)
    w = future.apply::future_apply(x.eval.grid, 1, FUN = weights_point,
                     x.design.grid = x.design.grid, h = h, m = m, del = del, ...,
                     future.seed = T)
    parallel::stopCluster(cl)
  }else{
    w = apply(x.eval.grid, 1, weights_point,
              x.design.grid = x.design.grid,
              h = h, m = m, del = del, ...)
  }

  if(del == 0){
    weights = w
  }else{
    if (grid.type == "less") {k = 2} else {k = 1}
    weights = slice_matrix(w, p*(p-1)/k, switch(del, 3, 6))
  }
  rm(w)
  L = list(design = x.design.grid, x.design = (1:p-1/2)/p, p = p, grid.type = grid.type,
           eval = x.eval.grid, x.eval = (1:p.eval-1/2)/p.eval, p.eval = p.eval, eval.type = eval.type,
           bandwidth = h, m = m, del = del, weights = weights)
}

#### Evaluation of the weights ####


#' Evaluates the weights from local_polynomial_weights with data
#'
#' @param W weights object from local_polynomial_weights
#' @param Z transformed data
#'
#' @return Matrix/Array with evaluation
#' @export
#'
#' @examples
#' 0 # TODO
eval_weights = function(W, Z){
  if(W$del == 0){
    if (W$grid.type == "less") {
      if (W$eval.type == "full") {
        M = matrix(0, W$p.eval, W$p.eval)
        M.up = upper.tri(M, T)
        M[M.up] = crossprod(W$weights, Z)
        M[!M.up] = t(M)[!M.up]
      } else {
        if(W$eval.type != "diagonal") {stop("eval.type needs to be diagonal or full")}
        M = crossprod(W$weights, Z)
      }
    } else {
      if(W$grid.type != "without diagonal") {stop("grid type not implemented")}
      if(W$p*(W$p-1) != length(Z)) {stop("Z is not in correct representation")}
      M = crossprod(W$weights, Z)
    }
    return(M)
  }else{
    if(W$eval.type == "full"){
      if(W$grid.type == "less"){
        weights_eval = apply(W$weights, 3, function(w){M = matrix(0, W$p.eval, W$p.eval)
        M.up = upper.tri(M, T)
        M[M.up] = crossprod(w, Z)
        M[!M.up] = t(M)[!M.up]
        M})
        weights_eval = array(weights_eval, dim = c(W$p.eval, W$p.eval, switch(W$del, 3, 6)))
      }else if(W$grid.type == "without diagonal"){
        weights_eval = apply(W$weights, 3, function(w){crossprod(w, Z)}) |>
          array(dim = c(W$p.eval, W$p.eval, switch(W$del, 3, 6)))
      }else{
        stop("grid type not implemented")
      }
    }else{
      if(W$eval.type != "diagonal"){stop("eval.type needs to be -diagonal- or -full-")}
      weights_eval = apply(W$weights, 3, function(w){crossprod(w, Z)})
    }
  }
  return(weights_eval)
}

