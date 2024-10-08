% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/local_polynomial_weights.R
\name{weights_point}
\alias{weights_point}
\title{Calculate the weights of bivariate local polynomial estimator at a single point}
\usage{
weights_point(x, x.design.grid, h, K = epak_2d, m = 1, del = 0)
}
\arguments{
\item{x}{design point as two dimensional vector at which the weights are calculated}

\item{x.design.grid}{grid of design points for calculation of the weights}

\item{h}{bandwidth parameter that needs to be chosen}

\item{K}{kernel function - default is two dimensional Epanechnikov kernel}

\item{m}{1 refers to local linear, 2 to local quadratic}

\item{del}{partial derivatives of which order shall be calculated?}
}
\value{
returns a matrix with one to six columns where the first column
 is the estimation of the actual function, the second is the partial derivative
 in the first direction (x), the third in direction 2 (y), the fourth xx, the
 fifth xy and the sixth in direction yy. In case `del = 0` its an row vector and other wise a matrix.
}
\description{
function that gets called from local.polynomial.weights
 to calculate the weights at a specific point. Therefore the function is not
 meant to be used singularly.
}
\examples{
weights_point(c(0.2, 0.3), observation_grid(15, comp = "less"), 0.2, del = 1)
}
