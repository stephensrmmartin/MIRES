##' Computes p_j - p_i not j, for each j.
##'
##' Assumes that mcmc contains params from ONE parameter across the indices whose differences are of interest (e.g, resid_random[1:K, 1]).
##' @title Outer subtraction for given params across MCMC samples.
##' @param mires Mires object.
##' @return MCMC Matrix of all differences. Posterior samples of all possible differences, minus duplicates.
##' @author Stephen R Martin
##' @keywords internal
.pairwise_diff_single <- function(mcmc) {
    stan_indices <- split_stannames(colnames(mcmc))$indices
    diffs <- apply(mcmc, 1, .sample_diff)
    colnames(diffs) <- .pairwise_diff_single_names(stan_indices[,1])
    diffs
}
##' @title Compute all differences of vector.
##' @param x Numeric vector (e.g., one row of mcmc samples)
##' @return Numeric Vector of lower-triangular outer difference matrix.
##' @author Stephen R Martin
##' @keywords internal
.sample_diff <- function(x) {
    out <- outer(x, x, "-")
    out[lower.tri(out)]
}
##' @title Generate labels for all differences of vector.
##' @param mcmcNames 
##' @return Character Vector of lower-triangular outer difference matrix. I.e., labels for .sample_diff.
##' @author Stephen R Martin
##' @keywords internal
.pairwise_diff_single_names <- function(mcmcNames) {
    out <- outer(mcmcNames, mcmcNames, function(x, y) {paste0(x, "-", y)})
    out <- out[lower.tri(out)]
    out
}
