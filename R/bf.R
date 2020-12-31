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

.density.logspline <- function(mcmc, lbound = 0, ...) {
    lso <- logspline(mcmc, lbound = lbound)
    df <- function(x) {
        dlogspline(x, lso)
    }
    return(df)
}

# TODO : Try this again, but with Weibull. Then, again with Stan/vb/optim.

.density.DP <- function(mcmc, iter = 500, mode = c("posterior", "est"), ...) {
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
predictMixture <- function(x, fit, K, pi = "pi", dens, params, R_params) {
    pi <- as.matrix(fit, pars = pi)
    params <- lapply(params, function(p){as.matrix(fit, pars = p)})
    names(params) <- R_params
    predfun <- function(x) {
        px <- rowSums(pi * matrix(do.call(dens, c(x = x, params)), nrow(pi), K))
        out <- c(mean = mean(px), sd = sd(px), quantile(px, c(.025, .975)))
        names(out)[3:4] <- c("Q2.5", "Q97.5")
        out
    }
    pxs <- cbind(x, t(sapply(x, predfun)))
    pxs
}
