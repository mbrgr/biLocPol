##### Mean functions for simulations #####
#' Mean Function
#'
#' @param x vector with points to evaluate function
#'
#' @return vector of same length as x
#' @export
#'
#' @examples
#' mu(0.3)
mu = function(x){
  x = 2*x - 1
  return(sin(3*pi*x)*exp(-2*abs(x)))
}



#' Generate FDA data
#' @description Generate full data from a process with Brownian motions and normal distributed errors
#'
#'
#' @param n number of observed curves
#' @param x.design vector with design points (length p)
#' @param f mean function
#' @param r.process random process
#' @param process.arg list with further arguments that can be shall bes passed to the process
#' @param r.error distribution of random errors
#' @param eps.arg  list with further arguments that can be shall bes passed to the errors
#' @param transformed shall the weights be transformed according to observation_transformation?
#' @param decomposed shall the decomposition into the parts be returned?
#' @param ... further arguments to be passed
#'
#' @return TODO
#' @importFrom stats rnorm
#' @export
#'
#' @examples
#' 0 # TODO
FDA_observation = function(n, x.design, f = mu,
                           r.process = BM, process.arg = list(),
                           r.error = rnorm, eps.arg = list(), transformed = F,
                           decomposed = F, ...){

  p = length(x.design)
  f.eval = f(x.design) # R^p --> not necessary, since it cancels out anyway
  process.arg = c(list(n = n, t = x.design), process.arg)
  process = do.call(r.process, process.arg) # R^{n x p}
  eps.arg = c(list(n = n*p), eps.arg)
  eps = matrix(do.call(r.error, eps.arg), n, p) # R^{n x p}
  Y = f.eval + process + eps
  rm(process.arg); rm(eps.arg)

  if(transformed)
    return(observation_transformation(Y, ...) )
  if(decomposed)
    return(list(f.eval, process, eps, Y))
  return(Y)
}


#' Brownian Motion
#'
#' @param n TODO
#' @param t TODO
#' @param sigma TODO
#'
#' @return TODO
#' @importFrom stats rnorm
#' @export
#'
#' @examples
#' 0 # TODO
BM = function(n, t, sigma = 1){
  v = diff(c(0, t))
  p = length(v)
  t(replicate(n, cumsum(rnorm(p, 0, sigma*v^0.5)) ))
}


#' Covariance kernel of Ornstein-Uhlenbeck process
#'
#' @param t 2 dimensional
#' @param sigma influence of brownian motion to process
#' @param theta mean reversion
#' @param sigma0 Variance of starting distribution, 0 if X0 is deterministic
#'
#' @return value of covariance
#' @export
#'
#' @examples
#' 0 # TODO
cov.ou = function(t, sigma, theta, sigma0 = 0){
  sigma^2/(2*theta) * (exp(-theta*abs(t[1]-t[2])) - exp(-theta*(t[1]+t[2]))) + sigma0^2*exp(-theta*(t[1]+t[2]))
}

#' Ornstein Uhlenbeck Process
#' @description Uses the ... package !!
#'
#' @param n how many processes
#' @param t points of evaluation
#' @param mu mean
#' @param alpha TODO
#' @param sigma TODO
#' @param x0 TODO
#'
#' @return TODO
#' @importFrom goffda r_ou
#' @export
#'
#' @examples
#' 0 # TODO
OU = function(n, t = seq(0, 1, len = 201), mu = 0, alpha = 1, sigma = 1, x0 = 0){
  goffda::r_ou(n, t, mu, alpha, sigma, x0)$data
}


#' Generates values to test functions
#'
#' @param ... grid type
#'
#' @return values
#' @export
#'
#' @examples
#' generate_test_data()
generate_test_data = function(...){
  n = 100
  p = 25
  x.design = (1:p - 0.5)/p
  Y = FDA_observation(n, x.design)
  list(
    n = n, p = p, x.design = x.design,
    p.eval = 50,
    h = 0.2,
    m = 1,
    del = 0,
    x = c(0.2, 0.3),
    Y = Y, Z = observation_transformation(Y, ...)
  )
}
