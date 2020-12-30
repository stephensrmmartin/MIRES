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
