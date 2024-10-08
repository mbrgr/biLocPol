% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data_generating_functions.R
\name{cov_z_2rv}
\alias{cov_z_2rv}
\title{cov_z_2rv}
\usage{
cov_z_2rv(t, a = 2/3, b = sqrt(2) * 2/3, c = 1.25)
}
\arguments{
\item{t}{Two dimensional evaluation vector of covariance kernel}

\item{a}{Parameter}

\item{b}{Parameter}

\item{c}{Parameter}
}
\value{
Single value of covariance kernel
}
\description{
cov_z_2rv
}
\examples{
cov_z_2rv(c(0.4, 0.3))
}
