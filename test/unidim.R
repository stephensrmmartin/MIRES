library(rstan)
library(MIRES)

set.seed(14)
## J <- 10
## K <- 20
## n <- 40
J <- 8
K <- 10
n <- 100

fixed <- list(lambda = rep(.7, J), resid_log = rep(log(sqrt(1 - .7^2)), J), nu = rep(0, J))
mipattern <- list("none")
etadist <- "std"

## mipattern <- list("none")

d <- MIRES:::datagen_uni(J, K, n, fixed, mipattern, etadist = "std")

ds <- d$df

fit_s_dep_marg <- mires(myLatent ~ x_1 + x_2 + x_3 + x_4 + x_5 + x_6 + x_7 + x_8, group = group, ds,
                       inclusion_model = "dep",
                       identification = "sum",
                       save_scores = FALSE,
                       prior = c(0, .25), iter = 1000)

fit_s_ind_marg <- mires(myLatent ~ x_1 + x_2 + x_3 + x_4 + x_5 + x_6 + x_7 + x_8, group = group, ds,
                       inclusion_model = "indep",
                       identification = "sum",
                       save_scores = FALSE,
                       prior = c(0, .25), iter = 1000)

fit_h_dep_marg <- mires(myLatent ~ x_1 + x_2 + x_3 + x_4 + x_5 + x_6 + x_7 + x_8, group = group, ds,
                       inclusion_model = "dep",
                       identification = "hier",
                       save_scores = FALSE,
                       prior = c(0, .25), iter = 1000)

fit_s_dep_cond <- mires(myLatent ~ x_1 + x_2 + x_3 + x_4 + x_5 + x_6 + x_7 + x_8, group = group, ds,
                       inclusion_model = "dep",
                       identification = "sum",
                       save_scores = TRUE,
                       prior = c(0, .25), iter = 1000)
