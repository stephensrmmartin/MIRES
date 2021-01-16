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
    funs <- apply(samps, 2, dlogspline)

    # Which failed?
    failed <- sapply(funs, function(x){!is.function(x)})
    failed <- which(failed)
    if(length(failed) > 0) {
        warning(length(failed), " logspline density approximations failed. Using (experimental) Dirichlet Process approximations for failed chains.")
    }

    # Recompute using HNormal DP
    funs[failed] <- sapply(samps[, failed], ddirichletprocess_spike, ...)

    return(funs)
}

savage_dickey <- function(posterior_dens_fun, prior_dens_fun, ...) {
    posterior_dens_fun(0) / do.call(prior_dens_fun, ...)
}
