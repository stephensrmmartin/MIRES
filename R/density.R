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
##' @keywords internal
rhmre <- function(n, mu = 0, sigma = 1) {
    # TODO : Extend this to allow the full HMRE structure (J items, P params, J*P uniques, 1 tau)
    tau <- rnorm(n, mu, sigma)
    item <- rnorm(n, mu, sigma)
    param <- rnorm(n, mu, sigma)
    unique <- rnorm(n, mu, sigma)

    hmre <- exp(tau + item + param + unique)

    abs(rnorm(n, 0, hmre))
}

phmre <- function(q, mu = 0, sigma = 1, lower.tail = TRUE) {
    f <- function(q, mu, sigma) {
        integrate(dhmre, 0, q, mu = mu, sigma = sigma)[["value"]]
    }
    f <- Vectorize(f)
    out <- f(q, mu, sigma)
    if(!lower.tail) out <- 1 - out
    out
}

##' @title Create logspline-based density function.
##' @param mcmc MCMC samples.
##' @param lbound Integer (Default: 0).
##' @param ... Not used.
##' @return Density Function.
##' @author Stephen Martin
##' @import logspline
##' @keywords internal
dlogspline <- function(mcmc, lbound = 0, ...) {
    lso <- logspline::logspline(mcmc, lbound = lbound, ...)
    df <- function(x) {
        logspline::dlogspline(x, lso)
    }
    if(class(lso) == "logspline") {
        return(df)
    } else {
        return(NA)
    }
}

##' @title Create dirichletprocess (exponential) based density function.
##' @param mcmc MCMC samples.
##' @param iter MH Iterations to run.
##' @param mode posterior or est.
##' @param ... Not used.
##' @return Function returning a matrix (if posterior) or vector (if est).
##' @author Stephen R. Martin
##' @import dirichletprocess
##' @keywords internal
ddirichletprocess <- function(mcmc, iter = 500, mode = c("posterior", "est"), ...) {
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
ddirichletprocess_stan <- function(mcmc, mode = "est", K = 200, model = "dpHNormal", ...) {
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
            predict_DP(x, stanOut, K,
                           dens = params[[model]]$dens,
                           params = params[[model]]$params,
                           R_params = params[[model]]$R_params)
        }
    } else if (mode == "est") {
        fun <- function(x) {
            predict_DP(x, stanOut, K,
                           dens = params[[model]]$dens,
                           params = params[[model]]$params,
                           R_params = params[[model]]$R_params)[,"mean"]
        }
    } 
     
    return(fun)
}

##' @title Create Stan-based spike-mixture DP based density estimation function.
##' @param mcmc MCMC samples.
##' @param mode posterior or est.
##' @param K Number of DP components (Default: 200)
##' @param spike_scale Numeric (Default: .00001). The scale of the half-normal spike.
##' @param ... Not used.
##' @return Function.
##' @author Stephen R Martin
##' @keywords internal
ddirichletprocess_spike <- function(mcmc, mode = "est", K = 200, spike_scale = .00001, ...) {
    dots <- list(...)
    stan_data <- list(N = length(mcmc),
                      y = mcmc,
                      K = K)
    stanOut <- rstan::vb(stanmodels[["dpHNormalSpike"]],
                         data = stan_data,
                         importance_resampling = FALSE,
                         tol_rel_obj = dots$tol_rel_obj %IfNull% .005)

    pi_mix <- as.matrix(stanOut, pars = c("pi_mix"))

    fun.samps <- function(x) {
        # Samples of the DP-part of prediction
        samps <- predict_DP(x, stanOut, K,
                        dens = dpnorm,
                        params = c("location", "scale"),
                        R_params = c("mu", "sigma"),
                        samps = TRUE)
        samps <- apply(samps, 2, function(x){x * (pi_mix)}) # Multiply DP_y by 1 - pi_mix
        spike_dens <- dpnorm(x, 0, spike_scale)
        spike_samps <- sapply(spike_dens, function(x) {x * (1 - pi_mix)})
        samps <- samps + spike_samps
        return(samps)
    }

    fun.posterior <- function(x) {
        samps <- fun.samps(x)
        out <- t(apply(samps, 2, function(x) {
            c(mean = mean(x), sd = sd(x), quantile(x, c(.025, .975)))
        }))
        colnames(out) <- c("mean", "sd", "Q2.5", "Q97.5")
        out
    }

    fun.est <- function(x) {
        out <- fun.posterior(x)[,"mean"]
        out
    }

    if(mode == "est") {
        return(fun.est)
    } else if(mode == "posterior") {
        return(fun.posterior)
    } else if(mode == "samps") {
        return(fun.samps)
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
##' @param param Data.frame of parameters corresponding to d and r distribution  functions. (E.g., data.frame(mean = rnorm(50), sd = abs(rnorm(50, 0, .5))))
##' @param alpha Numeric. The alpha parameter to the DP.
##' @param f Character. Root name of base or kernel function (e.g., "norm", "exp").
##' @return List of data (y), weights (pi), params (param), and the true density function (d).
##' @author Stephen Martin
##' @keywords internal
simulate_DP <- function(N, K, param, alpha, f) {
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
predict_DP <- function(x, fit, K, pi = "pi", dens, params, R_params, samps = FALSE) {
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
##' Computes the implied densities of random effect differences given HMRE prior.
##'
##' The HMRE prior for the RE-SD is \eqn{\int N^+(\sigma_p | exp(h_p))LN(h_p | 4\mu, \sqrt{4}\sigma)dh_p}.
##' The random effects are distributed as \eqn{u_{k,p} \sim N(0, \sigma_p)}.
##' The implied prior is therefore \eqn{u_{k,p} - u_{\lnot k, p} \sim N(0, \sqrt{2}\sigma)}.
##' Note that there is a singularity at 0, because the integrand at sigma = 0 is an infinite spike.
##' We currently integrate (using a change of variables) starting at machine precision-zero. Consider this the approximation of the limit as we approach 0 positively.
##' This is therefore divergent when assessed at a difference of zero, due to the RESD taking on a zero value (and an infinite function value).
##' This is expected, as the limit of a Gaussian as sigma -> 0 is the Dirac delta function.
##' @title Implied density for pairwise differences given HMRE prior.
##' @param x Numeric. Difference in random effects.
##' @param mu Numeric. HMRE Prior location.
##' @param sigma Numeric. (Default: 1; must be > 0). HMRE prior scale.
##' @return Numeric vector.
##' @author Stephen R. Martin
##' @keywords internal
dhmre_pairwise <- function(x, mu = 0, sigma = 1) {

    jointprior <- function(x, resd, mu, sigma) {
        dnorm(x, 0, sqrt(2) * resd) * dhmre(resd, mu, sigma)
    }
    g <- function(u) {
        0 + u / (1 - u)
    }
    g_jacobian <- function(u) {
        1 / ((1 - u)^2)
    }
    jointprior_CoV <- function(x, u, mu, sigma) {
        resd <- g(u)
        dnorm(x, 0, sqrt(2) * resd) * dhmre(resd, mu, sigma) * g_jacobian(u)
    }
    margprior <- function(x, mu, sigma) {
        #### Untransformed space; 
        ## pracma::quadinf(jointprior, xa = 0, xb = Inf, x = x, mu = mu, sigma = sigma)[["Q"]]
        ## integrate(jointprior, lower = 0, upper = Inf, x = x, mu = mu, sigma = sigma, stop.on.error=TRUE)[["value"]]
        #### Using a change of variables g(0) -> 0; g(1) -> Inf
        ## integrate(jointprior_CoV, lower = .Machine$double.eps, upper = 1, x = x, mu = mu, sigma = sigma, stop.on.error=FALSE, subdivisions = 1000)[["value"]]
        ## pracma::quadinf(jointprior_CoV, xa = 0, xb = 1, x = x, mu = mu, sigma = sigma)[["Q"]]
        cubature::cubintegrate(jointprior_CoV, .Machine$double.eps, 1, fDim = 1, method = "hcubature", x = x, mu = mu, sigma = sigma)[["integral"]]
    }
    margprior <- Vectorize(margprior)

    ## as.numeric(margprior(x, mu, sigma)["value",])
    margprior(x, mu, sigma)
    
}

rhmre_pairwise <- function(n, mu = 0, sigma = 1) {
    resd <- rhmre(n, mu, sigma)
    us <- t(sapply(resd, function(x){rnorm(2, 0, sqrt(2) * x)}))
    ds <- us[,1] - us[,2]
    ds
}
