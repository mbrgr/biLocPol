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
cov_ou = function(t, sigma, theta, sigma0 = 0){
  sigma^2/(2*theta) * (exp(-theta*abs(t[1]-t[2])) - exp(-theta*(t[1]+t[2]))) + sigma0^2*exp(-theta*(t[1]+t[2]))
}

#' Ornstein Uhlenbeck Process
#' @description Uses the [goffda] package.
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



#' z_2rv
#'
#' @param n Amount of observated curves in FDA setting
#' @param p Amount of Design points in FDA setting.
#' @param t Vector for design points if the design shall differ from the equidistant design. If not provided equidistant setting is used.
#' @param a Parameter
#' @param b Parameter
#' @param c Parameter
#'
#' @return Matrix of dimension n x p.
#' @export
#'
#' @examples z_2rv(n = 10)
z_2rv = function(n = 1, p = 100, t = NULL,
                 a = 2/3, b = sqrt(2)*2/3, c = 1.25){
  if( is.null(t) ) {
    t = (1:p - 1/2)/p
  }
  N = matrix(rnorm(2*n, 0, 1), n, 2)
  apply(N, 1, function(m){
    a * m[1] * sin(pi * t) +
      b * m[2] * cos(c*pi*t)
  }) |> t()
}


#' cov_z_2rv
#'
#' @param t Two dimensional evaluation vector of covariance kernel
#' @param a Parameter
#' @param b Parameter
#' @param c Parameter
#'
#' @return Single value of covariance kernel
#' @export
#'
#' @examples cov_z_2rv(c(0.4, 0.3))
cov_z_2rv = function(t, a = 2/3, b = sqrt(2)*2/3, c = 1.25){
  a^2*sin(pi*t[1])*sin(pi*t[2]) +
    b^2*cos(c*pi*t[1])*cos(c*pi*t[2])
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
