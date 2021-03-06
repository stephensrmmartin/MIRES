% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utility.R
\name{tidy_stanpars}
\alias{tidy_stanpars}
\title{Tidy up a vector of stan names into a data frame.}
\usage{
tidy_stanpars(stannames, labs = NULL, ...)
}
\arguments{
\item{stannames}{Character vector of stan names.}

\item{labs}{Optional. Character vector for what indices meaningfully pertain to. E.g., c("Factor", "Item") for an FxJ matrix.}

\item{...}{Optional. Named vectors of labels corresponding to indices. E.g., if labs = "Item", and you include \code{Item = colnames(model_matrix)}, then the indices (numeric) are replaced by the name in the supplied vector. E.g., the numeric i will be replaced by the i-th name.}
}
\value{
Data frame containing parameters and the (optionally named) indices.
}
\description{
Tidy up a vector of stan names into a data frame.
}
\author{
Stephen R Martin
}
\keyword{internal}
