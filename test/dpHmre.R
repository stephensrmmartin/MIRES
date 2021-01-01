library(rstan)
library(dirichletprocess)
library(logspline)

dposNorm <- function(x, mu, sd) {
    EnvStats::dnormTrunc(x, mu, sd, 0, Inf)
}

# Generate data
set.seed(13)
dpTest <- rhmre(2000, -1, 1)

# Compile stan model
dpWeibull <- stan_model("inst/stan/dpWeibull.stan")
dpExp <- stan_model("inst/stan/dpExp.stan")
dpHn <- stan_model("inst/stan/dpHNormal.stan")

# Estimate [sampling, vb, or optimizing]
K <- 50
stan_data <- list(N = length(dpTest), y = dpTest, K = K)
## dpOut_sample <- sampling(dpHn, data = stan_data, iter = 500, cores = 4)
dpOut <- vb(dpHn, data = stan_data,
            algorithm = "meanfield",
            importance_resampling = TRUE, tol_rel_obj = .005)

# Fit using R::dirichletprocess
## dpOut2 <- Fit(DirichletProcessWeibull(dpTest, c(3,2,2)), 500)
dpOut2 <- Fit(DirichletProcessExponential(dpTest), 500)

# Plot hist, true density, stan density, and dirichletprocess density
hist(dpTest, breaks = 300, probability = TRUE)
curve(dhmre(x,-1, 1), 0, 4, n = 1000, add = TRUE, lwd = 3)
curve(predictMixture(x, dpOut, K, dens = dposNorm, params = c("location", "scale"), R_params = c("mu", "sd"))[, "mean"], 0, 4, n = 1000, add = TRUE, col = "red", lwd = 3)
## curve(predictMixture(x, dpOut, K, dens = dexp, params = c("rate"), R_params = c("rate"))[,"mean"], 0, 4, n = 1000, add = TRUE, col = "red")
curve(dlogspline(x,logspline(dpTest, lbound = 0)), 0, 5, n = 1000, add = TRUE, col = "blue")
curve(PosteriorFunction(dpOut2)(x), 0, 4, n = 1000, add = TRUE, col = "green")
