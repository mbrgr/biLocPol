% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data_generating_functions.R
\name{FDA_observation}
\alias{FDA_observation}
\title{Generate FDA data}
\usage{
FDA_observation(
  n,
  x.design,
  f = mu,
  r.process = BM,
  process.arg = list(),
  r.error = rnorm,
  eps.arg = list(),
  transformed = F,
  decomposed = F,
  ...
)
}
\arguments{
\item{n}{number of observed curves}

\item{x.design}{vector with design points (length p)}

\item{f}{mean function}

\item{r.process}{random process}

\item{process.arg}{list with further arguments that can be shall bes passed to the process}

\item{r.error}{distribution of random errors}

\item{eps.arg}{list with further arguments that can be shall bes passed to the errors}

\item{transformed}{shall the weights be transformed according to observation_transformation?}

\item{decomposed}{shall the decomposition into the parts be returned?}

\item{...}{further arguments to be passed}
}
\value{
TODO
}
\description{
Generate full data from a process with Brownian motions and normal distributed errors
}
\examples{
0 # TODO
}
