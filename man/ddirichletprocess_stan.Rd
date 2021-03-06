% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/density.R
\name{ddirichletprocess_stan}
\alias{ddirichletprocess_stan}
\title{Create Stan-based density function.}
\usage{
ddirichletprocess_stan(mcmc, mode = "est", K = 200, model = "dpHNormal", ...)
}
\arguments{
\item{mcmc}{MCMC samples.}

\item{mode}{posterior or est.}

\item{K}{Number of DP components (Default: 200)}

\item{model}{dpHNormal, dpExp, dpGauss, or dpWeibull (Default: dpHNormal).}

\item{...}{Not used.}
}
\value{
Function returning vector (if est) or matrix (if posterior)
}
\description{
Create Stan-based density function.
}
\author{
Stephen Martin
}
\keyword{internal}
