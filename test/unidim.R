library(rstan)
library(MIRES)

set.seed(14)
J <- 8
K <- 5
n <- 50

fixed <- list(lambda = rep(.7, J), resid_log = rep(log(sqrt(1 - .7^2)), J), nu = rep(0, J))
etadist <- "std"

mipatterns <- list(
    none = list("none"),
    constant = list("constant", .3),
    items = list("items", .3, J/2),
    loadings = list("params", .3, 0),
    resid = list("params", .3, 1),
    nu = list("params", .3, 2)
)

datasets <- lapply(mipatterns, MIRES:::datagen_uni, J = J, K = K, n = n, etadist = etadist, fixed = fixed)

fits <- lapply(datasets, function(x) {
    df <- x$df
    fit <- mires(myFactor ~ x_1 + x_2 + x_3 + x_4 + x_5 + x_6 + x_7 + x_8, group, df, iter = 1000)
    fit
})

fits_summaries <- lapply(fits, function(x) {
    summary(x)
})

BF01s <- sapply(fits_summaries, function(x) {
    x$summary$resd$BF01
})

BFsls <- sapply(fits_summaries, function(x) {
    x$summary$resd$`BF(SD <= 0.1)`
})

round(BF01s, 3)
round(BFsls, 3)

fits_ind <- lapply(datasets, function(x) {
    df <- x$df
    fit <- mires(myFactor ~ x_1 + x_2 + x_3 + x_4 + x_5 + x_6 + x_7 + x_8, group, df, inclusion_model = "ind", iter = 1000)
    fit
})

fits_ind_summaries <- lapply(fits, summary)

BF01s_ind <- sapply(fits_ind_summaries, function(x) {
    x$summary$resd$BF01
})

BFsls_ind <- sapply(fits_ind_summaries, function(x) {
    x$summary$resd$`BF(SD <= 0.1)`
})
