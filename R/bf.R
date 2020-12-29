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
    tau <- rnorm(n, mu, sigma)
    item <- rnorm(n, mu, sigma)
    param <- rnorm(n, mu, sigma)
    unique <- rnorm(n, mu, sigma)

    hmre <- exp(tau + item + param + unique)

    abs(rnorm(n, 0, hmre))
}
