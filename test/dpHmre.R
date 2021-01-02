library(rstan)
library(dirichletprocess)
library(logspline)

dposNorm <- function(x, mu, sd) {
    EnvStats::dnormTrunc(x, mu, sd, 0, Inf)
}

# Generate data
set.seed(13)
dpTest <- MIRES:::rhmre(4000, -1, 1)

# Get Stan method
dfun.stan <- MIRES:::.density.stan(dpTest, K = 200)

# Compare
MIRES:::dhmre(0, -1, 1)
dfun.stan(0)
logspline::dlogspline(0, logspline::logspline(dpTest, lbound = 0))

# Plot hist, true density, stan density, and dirichletprocess density
hist(dpTest, breaks = 300, probability = TRUE)
curve(MIRES:::dhmre(x,-1, 1), 0, 4, n = 1000, add = TRUE, lwd = 3)
curve(dfun.stan(x), 0, 4, n = 1000, add = TRUE, col = "red", lwd = 3)
curve(dlogspline(x,logspline(dpTest, lbound = 0)), 0, 5, n = 1000, add = TRUE, col = "blue")
