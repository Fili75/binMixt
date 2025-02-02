% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mixt.R
\name{em_bin_gaus_comp}
\alias{em_bin_gaus_comp}
\title{Binned CL-EM algorithm for multivariate binned diagonal Gaussian mixture}
\usage{
em_bin_gaus_comp(bin, pi0, mu0, sigma0, tol, tim, print = F)
}
\arguments{
\item{bin}{Binned data}

\item{pi0}{Starting point for proportions}

\item{mu0}{Starting point for means}

\item{sigma0}{Starting point for variances}

\item{tol}{Tolerance of the algorithm}

\item{tim}{Maximum number of iterations}
}
\value{
List with original binned data, final estimates, log-likelihood and parameters sequences and medium time per iteration
}
\description{
Binned CL-EM algorithm for multivariate binned diagonal Gaussian mixture
}
