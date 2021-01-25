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

## fit_re <- mires(myfactor ~ x_1 + x_2 + x_3 + x_4 + x_5 + x_6 + x_7 + x_8 + x_9 + x_10, group = group, ds, iter = 1000, prior_only = FALSE, hmre_mu = 0, hmre_scale = .25, hmre = FALSE)
## fit_re_prior <- mires(myfactor ~ x_1 + x_2 + x_3 + x_4 + x_5 + x_6 + x_7 + x_8 + x_9 + x_10, group = group, ds, iter = 1000, prior_only = TRUE, hmre_mu = 0, hmre_scale = .25, hmre = FALSE)
fit_re <- readRDS("~/Output/MIRES/fit_re.Rds")
## saveRDS(fit_re, "~/Output/MIRES/fit_re.Rds")

# There doesn't seem to be any big difference between these; wtf. Maybe due to the data?
# No difference really when "items", .4, 5/10, K = 20, J = 10, n = 40

# Try again with different settings

set.seed(14)
J <- 10
K <- 5
n <- 100

fixed <- list(lambda = rep(.7, J), resid_log = rep(log(sqrt(1 - .7^2)), J), nu = rep(0, J))
mipattern <- list("none")
etadist <- "std"

d <- MIRES:::datagen_uni(J, K, n, fixed, mipattern, etadist = etadist)
ds <- d$df


## fit <- mires(myfactor ~ x_1 + x_2 + x_3 + x_4 + x_5 + x_6 + x_7 + x_8 + x_9 + x_10, group = group, ds, iter = 1000, prior_only = FALSE, hmre_mu = 0, hmre_scale = .25)
## fit_re <- mires(myfactor ~ x_1 + x_2 + x_3 + x_4 + x_5 + x_6 + x_7 + x_8 + x_9 + x_10, group = group, ds, iter = 1000, prior_only = FALSE, hmre_mu = 0, hmre_scale = .25, hmre = FALSE)
## saveRDS(fit, "~/Output/MIRES/fit_small.Rds")
## saveRDS(fit_re, "~/Output/MIRES/fit_re_small.Rds")
fit <- readRDS("~/Output/MIRES/fit_small.Rds")
fit_re <- readRDS("~/Output/MIRES/fit_re_small.Rds")
