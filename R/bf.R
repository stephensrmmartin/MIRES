##' For each RE-SD, approximates the marginal posterior density from MCMC samples for use in BF calculations.
##'
##' Starts by computing (lower-bounded) logspline approximations.
##' If these fail, it uses the Dirichlet process with positive-normal kernels as an approximation.
##' @title Create marginal posterior density function approximations for random effect SDs
##' @param mires mires object.
##' @param add_zero Logical (Default: TRUE). Whether to add a zero to samples.
##' @param ... Args passed onto .density.stan.
##' @return List of approximate density functions.
##' @author Stephen Martin
##' @keywords internal
posterior_density_funs_sigmas <- function(mires, add_zero = TRUE, ...) {
    pars <- "random_sigma"
    samps <- as.matrix(mires$fit, pars = pars)
    n_cols <- ncol(samps)
    if(add_zero) {
        samps <- rbind(0, samps)
    }

    # Try logspline FIRST, and fix with DP if failed.
    funs <- apply(samps, 2, dlogspline, error.action = 1)

    # Which failed?
    failed <- sapply(funs, function(x){!is.function(x)})
    failed <- which(failed)
    if(length(failed) > 0) {
        warning(length(failed), " logspline density approximations failed. Using (experimental) Dirichlet Process approximations for failed chains.")
    }

    # Recompute using HNormal DP
    funs[failed] <- sapply(as.data.frame(samps[, failed]), ddirichletprocess_spike, ...)

    return(funs)
}

savage_dickey <- function(posterior_dens_fun, prior_dens_fun, ...) {
    posterior_dens_fun(0) / prior_dens_fun(...)
}

prob_less_than <- function(mcmc, less_than) {
    mean(mcmc <= less_than)
}
##' Computes the BF12, where 1 is less than and 2 is greater than.
##'
##' The BF12 here is BF for a parameter being less than a threshold, t, vs the parameter being greater than t.
##' This borrows from the encompassing approach, where u is the unconstrained prior:
##' BF1u = p(D|H = 1) / p(D|H = u)
##' BF2u = p(D|H = 2) / p(D|H = u)
##' BF12 = BF1u / BF2u.
##'
##' BF1u = int p(D|H = 1, theta_1)p(theta_1 | H = 1) / p(D|H = u, theta_u)p(theta_u | H=u)
##' 
##' @title Compute BF(Less than)
##' @param mcmc MCMC vector.
##' @param less_than Value to test.
##' @param prior_cumul_fun CDF function.
##' @param ... 
##' @return BF12 value.
##' @author Stephen Martin
##' @keywords internal
bflt <- function(mcmc, less_than, prior_cumul_fun, ...) {
    post_prob_lt <- prob_less_than(mcmc, less_than)
    prior_prob_lt <- prior_cumul_fun(less_than, ...)
    bf1u <- post_prob_lt / prior_prob_lt

    post_prob_gt <- 1 - post_prob_lt
    prior_prob_gt <- 1- prior_prob_lt
    bf2u <- post_prob_gt / prior_prob_gt
    bf12_out <- bf1u / bf2u

    bf12_out
}
