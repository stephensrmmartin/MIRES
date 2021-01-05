##' @title Density for hmre prior on RE SDs.
##' @param x Numeric.
##' @param mu Numeric. HMRE Prior location.
##' @param sigma Numeric. (Default: 1; must be > 0). HMRE prior scale.
##' @return Density.
##' @author Stephen R. Martin
##' @keywords internal
dhmre <- function(x, mu = 0, sigma = 1) {
    # SD[p] ~ N^+ (0, exp(tau + item(p) + param(p) + unique(p)))
    # tau, item(p), param(p), unique(p) ~ N(mu, sigma)
    # tau + item(p) + param(p) + unique(p) ~implied~ N(0, 2*sigma)
    LN_sigma <- 2 * sigma
    LN_mu <- 4 * mu
    hn_given_hmre <- function(x, LN_mu, LN_sigma, hmre) {
        dnorm(x, 0, hmre) * 2 * dlnorm(hmre, LN_mu, LN_sigma) # Half normal
    }
    int_hn_over_hmre <- function(x, LN_mu, LN_sigma) {
        integrate(hn_given_hmre, lower = 0, upper = Inf, LN_mu = LN_mu, LN_sigma = LN_sigma, x = x)
    }
    int_hn_over_hmre <- Vectorize(int_hn_over_hmre)

    densities <- int_hn_over_hmre(x, LN_mu, LN_sigma)
    return(as.numeric(densities["value",]))
}

##' @title Random sampling from hmre prior on RE SDs.
##' @param n Integer.
##' @param mu Numeric. HMRE Prior location.
##' @param sigma Numeric. (Default: 1; must be > 0). HMRE prior scale.
##' @return Vector.
##' @author Stephen R. Martin
rhmre <- function(n, mu = 0, sigma = 1) {
    # TODO : Extend this to allow the full HMRE structure (J items, P params, J*P uniques, 1 tau)
    tau <- rnorm(n, mu, sigma)
    item <- rnorm(n, mu, sigma)
    param <- rnorm(n, mu, sigma)
    unique <- rnorm(n, mu, sigma)

    hmre <- exp(tau + item + param + unique)

    abs(rnorm(n, 0, hmre))
}

# TODO : HDInterval method [and add a 0 to each MCMC sample (or just RE-SDs)]

# TODO : Add logspline method (GP-based density estimator? Or Dirichlet-process density estimator?)

# TODO : Add pairwise-hmre density estimator; implied prior over hmre prior of u[k] - u[not_k]
##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##' @title Create logspline-based density function.
##' @param mcmc MCMC samples.
##' @param lbound Integer (Default: 0).
##' @param ... Not used.
##' @return Density Function.
##' @author Stephen Martin
##' @keywords internal
.density.logspline <- function(mcmc, lbound = 0, ...) {
    lso <- logspline(mcmc, lbound = lbound)
    df <- function(x) {
        dlogspline(x, lso)
    }
    return(df)
}

# TODO : Try this again, but with Weibull. Then, again with Stan/vb/optim.
##' @title Create dirichletprocess (exponential) based density function.
##' @param mcmc MCMC samples.
##' @param iter MH Iterations to run.
##' @param mode posterior or est.
##' @param ... Not used.
##' @return Function returning a matrix (if posterior) or vector (if est).
##' @author Stephen R. Martin
##' @keywords internal
.density.dirichletprocess <- function(mcmc, iter = 500, mode = c("posterior", "est"), ...) {
    dpo <- dirichletprocess::DirichletProcessExponential(mcmc, ...)
    dpo <- dirichletprocess::Fit(dpo, iter)
    if(mode == "est") {
        return(dirichletprocess::PosteriorFunction(dpo))
    } else if (mode == "posterior") {
        return(function(x, prob = .95) {
            dirichletprocess::PosteriorFrame(dpo, x, ci_size = 1 - prob, ndraws = iter)
        })
    }
}

##' @title Create Stan-based density function.
##' @param mcmc MCMC samples.
##' @param mode posterior or est.
##' @param K Number of DP components (Default: 200)
##' @param model dpHNormal, dpExp, dpGauss, or dpWeibull (Default: dpHNormal).
##' @param ... Not used.
##' @return Function returning vector (if est) or matrix (if posterior)
##' @author Stephen Martin
##' @keywords internal
.density.stan <- function(mcmc, mode = "est", K = 200, model = "dpHNormal", ...) {
    dots <- list(...)
    stan_data <- list(N = length(mcmc),
                      y = mcmc,
                      K = K)
    stanOut <- rstan::vb(stanmodels[[model]],
                         data = stan_data,
                         importance_resampling = TRUE,
                         tol_rel_obj = dots$tol_rel_obj %IfNull% .005
                         )
    params <- list(
        dpHNormal = list(params = c("location", "scale"),
                         R_params = c("mu", "sigma"),
                         dens = dpnorm),
        dpExp = list(params = c("rate"),
                     R_params = "rate",
                     dens = dexp),
        dpGauss = list(params = c("mu", "sigma"),
                       R_params = c("mean", "sd"),
                       dens = dnorm),
        dpWeibull = list(params = c("shape", "scale"),
                         R_params = c("shape", "scale"),
                         dens = dweibull)
    )
    if(mode == "posterior") {
        fun <- function(x) {
            predictMixture(x, stanOut, K,
                           dens = params[[model]]$dens,
                           params = params[[model]]$params,
                           R_params = params[[model]]$R_params)
        }
    } else if (mode == "est") {
        fun <- function(x) {
            predictMixture(x, stanOut, K,
                           dens = params[[model]]$dens,
                           params = params[[model]]$params,
                           R_params = params[[model]]$R_params)[,"mean"]
        }
    } 
     
    return(fun)
}

.density.stan_spike <- function(mcmc, mode = "est", K = 200, spike_scale = .00001, ...) {
    dots <- list(...)
    stan_data <- list(N = length(mcmc),
                      y = mcmc,
                      K = K)
    stanOut <- rstan::vb(stanmodels[["dpHNormalSpike"]],
                         data = stan_data,
                         importance_resampling = TRUE,
                         tol_rel_obj = dots$tol_rel_obj %IfNull% .005)

    pi_mix <- as.matrix(stanOut, pars = c("pi_mix"))

    fun <- function(x) {
        # Samples of the DP-part of prediction
        samps <- predictMixture(x, stanOut, K,
                        dens = dpnorm,
                        params = c("location", "scale"),
                        R_params = c("mu", "sigma"),
                        samps = TRUE)
    }
}

dpnorm <- function(x, mu = 0, sigma = 1) {
    truncnorm::dtruncnorm(x, mean = mu, sd = sigma, a = 0, b = Inf)
}

rpnorm <- function(n, mu = 0, sigma = 1) {
    truncnorm::rtruncnorm(n, mean = mu, sd = sigma, a = 0, b = Inf)
}

##' @title Generate Truncated Dirichlet Process Mixture.
##' @param N Number of data points.
##' @param K Max cluster.
##' @param param Data.frame of parameters corresponding to d and r distribution  functions.
##' @param alpha Numeric. The alpha parameter to the DP.
##' @param f Character. Root name of base or kernel function (e.g., "norm", "exp").
##' @return List of data (y), weights (pi), params (param), and the true density function (d).
##' @author Stephen Martin
##' @keywords internal
##' @examples
##' y <- genMixture(1000, 50, data.frame(mean = rnorm(50), sd = abs(rnorm(50, 0, .5))), .4, "norm")
genMixture <- function(N, K, param, alpha, f) {
    pi <- genStickBreakPi(K, alpha)

    y <- numeric(N)
    dens <- eval(as.symbol(paste0("d", f)))
    rng <- eval(as.symbol(paste0("r",f)))

    k_i <- sample(1:K, N, TRUE, pi)
    ## y <- mapply(rng, param[k_i, ])
    y <- do.call(rng, c(n = N, param[k_i,, drop = FALSE]))

    densfun <- function(y) {
        sum(pi * do.call(dens, c(x = y, param)))
    }
    densfun <- Vectorize(densfun, "y")

    list (y = y, pi = pi, param = param, d = densfun)

}

##' @title Stick-breaking function.
##' @param K Max cluster.
##' @param alpha Numeric. The alpha parameter to the DP.
##' @return Vector of DP weights.
##' @author Stephen Martin
##' @keywords internal
genStickBreakPi <- function(K, alpha) {
    stick_slices <- rbeta(K - 1, 1, alpha)
    pi <- numeric(length = K)
    pi[1] <- stick_slices[1]
    for(k in 2:(K - 1)) {
        pi[k] <- stick_slices[k] * prod(1 - stick_slices[1:(k - 1)])
    }
    pi[k] <- prod(1 - stick_slices[1:(K - 1)])
    return(pi)
}
##' @title Prediction for DP density estimation models.
##' @param x Values for prediction.
##' @param fit Stan DP fit.
##' @param K Max cluster.
##' @param pi Character. Name of stan variable corresponding to DP weights.
##' @param dens Function. The density base function used.
##' @param params Character vector. Names of base/kernel function parameters in Stan (e.g., mu, sigma for normal base functions).
##' @param R_params Character vector. Names of corresponding parameters for the R equivalent (e.g., mean, sd in dnorm).
##' @return Matrix of posterior mean, sd, .025, and .975 intervals.
##' @author Stephen R. Martin
##' @keywords internal
predictMixture <- function(x, fit, K, pi = "pi", dens, params, R_params, samps = FALSE) {
    pi <- as.matrix(fit, pars = pi)
    params <- lapply(params, function(p){as.matrix(fit, pars = p)})
    names(params) <- R_params
    if(!samps) {
        predfun <- function(x) {
            px <- rowSums(pi * matrix(do.call(dens, c(x = x, params)), nrow(pi), K))
            out <- c(mean = mean(px), sd = sd(px), quantile(px, c(.025, .975)))
            names(out)[3:4] <- c("Q2.5", "Q97.5")
            out
        }
        pxs <- cbind(x, t(sapply(x, predfun)))
        return(pxs)
    } else if(samps) {
        predfun <- function(x) {
            px <- rowSums(pi * matrix(do.call(dens, c(x = x, params)), nrow(pi), K))
            px
        }
        pxs <- sapply(x, predfun)
        return(pxs)
    }
}

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
    funs <- apply(samps, 2, .density.logspline)

    # Which failed?
    failed <- which(!is.function(funs))
    if(length(failed) > 0) {
        warning(length(failed), " logspline density approximations failed. Using (experimental) Dirichlet Process approximations for failed chains.")
    }

    # Recompute using HNormal DP
    funs[failed] <- apply(samps[, failed], .density.stan, ...)

    return(funs)
}
