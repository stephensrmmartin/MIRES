% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/MIRES-package.R
\docType{package}
\name{MIRES-package}
\alias{MIRES-package}
\alias{MIRES}
\title{The 'MIRES' package.}
\description{
Estimates random effect latent measurement models, wherein the loadings, residual variances, intercepts, latent means, and latent variances all vary across groups. The random effect variances of the measurement parameters are then modeled using a hierarchical inclusion model, wherein the inclusion of the variances (i.e., whether it is effectively zero or non-zero) is informed by similar parameters (of the same type, or of the same item). This additional hierarchical structure allows the evidence in favor of partial invariance to accumulate more quickly, and yields more certain decisions about measurement invariance.
}
\references{
Stan Development Team (2020). RStan: the R interface to Stan. R package version 2.21.1. https://mc-stan.org

Martin, S. R., Williams, D. R., & Rast, P. (2019, June 18). Measurement Invariance Assessment with Bayesian Hierarchical Inclusion Modeling. <doi:10.31234/osf.io/qbdjt>
}
\seealso{
Useful links:
\itemize{
  \item Report bugs at \url{https://github.com/stephenSRMMartin/MIRES/issues}
}

}
\author{
\strong{Maintainer}: Stephen Martin \email{stephenSRMMartin@gmail.com} (\href{https://orcid.org/0000-0001-8085-2390}{ORCID})

Authors:
\itemize{
  \item Philippe Rast \email{rast.ph@gmail.com} (\href{https://orcid.org/0000-0003-3630-6629}{ORCID})
}

}
