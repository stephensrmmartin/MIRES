library(rstan)
library(MIRES)

set.seed(13)
J <- 10
K <- 20
n <- 40

fixed <- list(lambda = rep(.7, J), resid_log = rep(log(sqrt(1 - .7^2)), J), nu = rep(0, J))
mipattern <- list("items", .4, floor(J/2))
etadist <- "std"

mipattern <- list("none")

d <- MIRES:::datagen_uni(J, K, n, fixed, mipattern, etadist = NULL)

ds <- d$df

fit <- mires(myfactor ~ x_1 + x_2 + x_3 + x_4 + x_5 + x_6 + x_7 + x_8 + x_9 + x_10, group = group, ds, iter = 1000, prior_only = FALSE)

samps <- as.matrix(fit$fit, pars = "random_sigma")
dfuns <- apply(samps, 2, MIRES:::dlogspline)
dfuns_spike <- apply(samps, 2, MIRES:::ddirichletprocess_spike, K = 100, tol_rel_obj = .01)
dfuns_HN <- apply(samps, 2, MIRES:::ddirichletprocess_stan, K = 100, tol_rel_obj = .01)
