library(rstan)
library(dirichletprocess)

genMixtureGaussian <- function(N, K, mu, sigma, alpha) {
    pi <- genStickBreakPi(K, alpha)

    y <- numeric(length = N)
    ## for(i in 1:N) {
    ##     y[i] <- sum(pi * rnorm(K, mu, sigma))
    ## }
    for(i in 1:N) {
        k_i <- sample(1:K, 1, prob = pi)
        y[i] <- rnorm(1, mu[k_i], sigma[k_i])
    }

    densfun <- function(y) {
        sum(dnorm(y, mu, sigma) * pi)
    }
    densfun <- Vectorize(densfun, "y")
    list(y = y, pi = pi, mu = mu, sigma = sigma, d = densfun)
}

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

predictMixtureGaussian <- function(x, fit, K) {
    samps <- as.matrix(fit, pars = c("pi","mu","sigma"))
    predfun <- function(x) {
        pis <- samps[, paste0("pi[",1:K,"]")]
        mus <- samps[, paste0("mu[",1:K,"]")]
        sigmas <- samps[, paste0("sigma[",1:K,"]")]
        px <- rowSums(matrix(dnorm(x, mus, sigmas), nrow(pis), K) * pis)
        c(mean = mean(px), sd = sd(px), Q2.5 = quantile(px, .025), Q97.5 = quantile(px,.975))
    }
    pxs <- cbind(x, t(sapply(x, predfun)))
    colnames(pxs)[4:5] <- c("Q2.5","Q97.5")
    pxs
}


# Generate data
set.seed(14)
dpTest <- genMixtureGaussian(1000, 50, rnorm(50,0,3), abs(rnorm(50,0,.3)), 1)

# Compile stan model
dpGauss <- stan_model("test/dpGauss.stan")

# Estimate [sampling, vb, or optimizing]
stan_data <- list(N = length(dpTest$y), y = dpTest$y, K = 200)
## dpOut <- sampling(dpGauss, data = stan_data, cores = 4, iter = 1000)
dpOut <- vb(dpGauss, data = stan_data, algorithm = "meanfield", importance_resampling = TRUE, tol_rel_obj = .002)
## dpOut <- optimizing(dpGauss, data = stan_data, algorithm = "LBFGS", draws = 1000, hessian = TRUE)

# Fit using R::dirichletprocess
dpOut2 <- Fit(DirichletProcessGaussian(dpTest$y), 1000)

# Plot hist, true density, stan density, and dirichletprocess density
hist(dpTest$y, breaks = 100, probability = TRUE)
curve(dpTest$d(x), -10, 5, n = 1000, add = TRUE)
curve(predictMixtureGaussian(x,dpOut, 200)[,"mean"], -10, 5, n = 1000, add = TRUE, col = "red")
curve(predictMixtureGaussian(x,dpOut, 200)[,"Q2.5"], -10, 5, n = 1000, add = TRUE, col = "red", lty = "dashed")
curve(predictMixtureGaussian(x,dpOut, 200)[,"Q97.5"], -10, 5, n = 1000, add = TRUE, col = "red", lty = "dashed")
curve(PosteriorFunction(dpOut2)(x), -10, 5, n = 1000, add = TRUE, col = "green")
