% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mixt.R
\name{binmixt.compinit}
\alias{binmixt.compinit}
\title{CL-EM with Initialization axes by axes}
\usage{
binmixt.compinit(
  data,
  class,
  R,
  maxit,
  tol,
  nrep,
  maxit1,
  print = T,
  only.init = F
)
}
\arguments{
\item{data}{Needs to a matrix of data}

\item{class}{Number of groups of the mixture models}

\item{R}{Needs to be a matrix or a numerical vector (vector at moment only option). The element (i,j) indicates the density of the j-th grid along the i-th dimension}

\item{maxit}{Maximum number of iterations small-EM}

\item{tol}{Tolerance of EM}

\item{nrep}{Number of different initialization of small-EM}

\item{maxit1}{Maximum number of iterations final EM}
}
\value{
Binned CL-EM algorithm output
}
\description{
CL-EM with Initialization axes by axes
}
