% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mixt.R
\name{binmixt}
\alias{binmixt}
\title{Binned CL-EM algorithm}
\usage{
binmixt(data, class, R, maxit, tol, nrep, print = T)
}
\arguments{
\item{data}{Needs to a matrix of data}

\item{class}{Number of groups of the mixture models}

\item{R}{Needs to be a matrix or a numerical vector. The element (i,j) indicates the density of the j-th grid along the i-th dimension}

\item{maxit}{Maximum number of iterations EM}

\item{tol}{Tolerance of EM}

\item{nrep}{Number of different initialization of EM}
}
\value{
Binned CL-EM algorithm output
}
\description{
Binned CL-EM algorithm
}
