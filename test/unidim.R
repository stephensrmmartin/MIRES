library(rstan)
library(MIRES)

set.seed(14)
## J <- 10
## K <- 20
## n <- 40
J <- 10
K <- 3
n <- 100

fixed <- list(lambda = rep(.7, J), resid_log = rep(log(sqrt(1 - .7^2)), J), nu = rep(0, J))
mipattern <- list("items", .4, floor(J/2))
mipattern <- list("params", .3, 2)
etadist <- "std"

## mipattern <- list("none")

d <- MIRES:::datagen_uni(J, K, n, fixed, mipattern, etadist = NULL)

ds <- d$df

fit_hier_incl <- mires(myLatent ~ x_1 + x_2 + x_3 + x_4 + x_5 + x_6 + x_7 + x_8 + x_9 + x_10, group = group, ds,
                       inclusion_model = "dep",
                       identification = "hier",
                       save_scores = FALSE,
                       prior = c(0, .25), iter = 1000)


fit_sum_incl <- mires(myLatent ~ x_1 + x_2 + x_3 + x_4 + x_5 + x_6 + x_7 + x_8 + x_9 + x_10, group = group, ds,
                       inclusion_model = "dep",
                       identification = "sum",
                       save_scores = FALSE,
                       prior = c(0, .25), iter = 1000)
