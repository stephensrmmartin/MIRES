library(rstan)
library(MIRES)

set.seed(14)
J <- 10
K <- 20
n <- 40

fixed <- list(lambda = rep(.7, J), resid_log = rep(log(sqrt(1 - .7^2)), J), nu = rep(0, J))
mipattern <- list("items", .4, floor(J/2))
etadist <- "std"

## mipattern <- list("none")

d <- readRDS("~/Output/MIRES/fit_gen.Rds")
## d <- MIRES:::datagen_uni(J, K, n, fixed, mipattern, etadist = etadist)
## saveRDS(d, "~/Output/MIRES/fit_gen.Rds")

ds <- d$df

## fit <- mires(myfactor ~ x_1 + x_2 + x_3 + x_4 + x_5 + x_6 + x_7 + x_8 + x_9 + x_10, group = group, ds, iter = 1000, prior_only = FALSE)
## fit <- mires(myfactor ~ x_1 + x_2 + x_3 + x_4 + x_5 + x_6 + x_7 + x_8 + x_9 + x_10, group = group, ds, iter = 1000, prior_only = FALSE, hmre_mu = 0, hmre_scale = .25)
## fit_prior <- mires(myfactor ~ x_1 + x_2 + x_3 + x_4 + x_5 + x_6 + x_7 + x_8 + x_9 + x_10, group = group, ds, iter = 1000, prior_only = TRUE, hmre_mu = 0, hmre_scale = .25)
## saveRDS(fit, "~/Output/MIRES/fit.Rds")
fit <- readRDS("~/Output/MIRES/fit.Rds")

samps <- as.matrix(fit$fit, pars = "random_sigma")
samps[1,] <- 0 # add zero

# Try posterior_sd fn
funs <- MIRES:::posterior_density_funs_sigmas(fit)
sapply(funs, function(x){x(0)})

column <- 26

fun.hn <- MIRES:::ddirichletprocess_spike(samps[,column], K = 100, tol_rel_obj = .001)
fun.exp <- MIRES:::ddirichletprocess_stan(samps[,column], K = 100, tol_rel_obj = .001, model = "dpExp")
fun.wei <- MIRES:::ddirichletprocess_stan(samps[,column], K = 100, tol_rel_obj = .001, model = "dpWeibull")
fun.hn2 <- MIRES:::ddirichletprocess_stan(samps[,column], K = 100, tol_rel_obj = .001)

hist(samps[,column], probability = TRUE, breaks = 200)
curve(funs[[column]](x), add = TRUE, col = "red", n = 1000, 0, .2)
curve(fun.hn, add = TRUE, col = "green", n = 1000, 0, .2)
curve(fun.exp, add = TRUE, col = "blue",lty = "dashed", n = 1000, 0, .2)
curve(fun.wei, add = TRUE, col = "blue",lty = "dotted", n = 1000, 0, .2)
curve(fun.hn2, add = TRUE, col = "blue", n = 1000, 0, .2)

