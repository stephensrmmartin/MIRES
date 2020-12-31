library(rstan)
library(dirichletprocess)
library(logspline)

# Generate data
set.seed(15)
dpTest <- genMixture(1000, 50, data.frame(mean = rnorm(50,0,3), sd = abs(rnorm(50,0,.2))), 1,"norm")

# Compile stan model
dpGauss <- stan_model("test/dpGauss.stan")

# Estimate [sampling, vb, or optimizing]
stan_data <- list(N = length(dpTest$y), y = dpTest$y, K = 200)
dpOut <- vb(dpGauss, data = stan_data, algorithm = "meanfield", importance_resampling = TRUE, tol_rel_obj = .005)

# Fit using R::dirichletprocess
dpOut2 <- Fit(DirichletProcessGaussian(dpTest$y), 1000)

# Plot hist, true density, stan density, and dirichletprocess density
hist(dpTest$y, breaks = 100, probability = TRUE)
curve(dpTest$d(x), -10, 5, n = 1000, add = TRUE)
curve(predictMixture(x,dpOut, 200, dens = dnorm, params = c("mu","sigma"), R_params = c("mean","sd"))[,"mean"], -10, 5, n = 1000, add = TRUE, col = "red")
curve(PosteriorFunction(dpOut2)(x), -10, 5, n = 1000, add = TRUE, col = "green")
curve(dlogspline(x,logspline(dpTest$y)), -10, 5, n = 1000, add = TRUE, col = "purple")
