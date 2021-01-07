library(rstan)
library(dirichletprocess)
library(logspline)

dposNorm <- function(x, mu, sd) {
    EnvStats::dnormTrunc(x, mu, sd, 0, Inf)
}

# Generate data
set.seed(13)
hmre_mu <- -2
hmre_sigma <- .5
dpTest <- MIRES:::rhmre(4000, hmre_mu, hmre_sigma)

# Get Stan method
dfun.stan <- MIRES:::.density.stan(dpTest, K = 100)
dfun.stan_spike <- MIRES:::.density.stan_spike(dpTest, K = 100)

# Get spiked variant
spiked <- stan_model("inst/stan/dpHNormalSpike.stan")
stan_data <- list(N = length(dpTest), y = dpTest, K = 100)
spikedOut <- vb(spiked, stan_data, tol_rel_obj = .005)

# Compare
MIRES:::dhmre(0, hmre_mu, hmre_sigma)
dfun.stan(0)
dfun.stan_spike(0)
summary(spikedOut, pars = c("py_0", "pi_mix"))$summary
logspline::dlogspline(0, logspline::logspline(dpTest, lbound = 0))

# Plot hist, true density, stan density, and dirichletprocess density
hist(dpTest, breaks = 300, probability = TRUE)
curve(MIRES:::dhmre(x,hmre_mu, hmre_sigma), 0, .01, n = 1000, add = TRUE, lwd = 3)
curve(dfun.stan(x), 0, .01, n = 1000, add = TRUE, col = "red", lwd = 3)
curve(dfun.stan_spiked(x), 0, .01, n = 1000, add = TRUE, col = "red", lwd = 3)
curve(dlogspline(x,logspline(dpTest, lbound = 0)), 0, .01, n = 1000, add = TRUE, col = "blue")
