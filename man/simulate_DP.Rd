% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/density.R
\name{simulate_DP}
\alias{simulate_DP}
\title{Generate Truncated Dirichlet Process Mixture.}
\usage{
simulate_DP(N, K, param, alpha, f)
}
\arguments{
\item{N}{Number of data points.}

\item{K}{Max cluster.}

\item{param}{Data.frame of parameters corresponding to d and r distribution  functions. (E.g., data.frame(mean = rnorm(50), sd = abs(rnorm(50, 0, .5))))}

\item{alpha}{Numeric. The alpha parameter to the DP.}

\item{f}{Character. Root name of base or kernel function (e.g., "norm", "exp").}
}
\value{
List of data (y), weights (pi), params (param), and the true density function (d).
}
\description{
Generate Truncated Dirichlet Process Mixture.
}
\author{
Stephen Martin
}
\keyword{internal}
