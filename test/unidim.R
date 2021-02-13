library(rstan)
library(MIRES)

set.seed(14)
J <- 8
K <- 5
n <- 50

fixed <- list(lambda = rep(.7, J), resid_log = rep(log(sqrt(1 - .7^2)), J), nu = rep(0, J))
etadist <- NULL

mipatterns <- list(
    none = list("none"),
    constant = list("constant", .4),
    items = list("items", .4, J/2),
    loadings = list("params", .4, 0),
    resid = list("params", .4, 1),
    nu = list("params", .4, 2)
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

fits_ind_summaries <- lapply(fits_ind, summary)

BF01s_ind <- sapply(fits_ind_summaries, function(x) {
    x$summary$resd$BF01
})

BFsls_ind <- sapply(fits_ind_summaries, function(x) {
    x$summary$resd$`BF(SD <= 0.1)`
})

round(BF01s_ind, 3)
round(BFsls_ind, 3)

plot(log(BF01s), log(BF01s_ind))
abline(0,1)

plot(log(BFsls), log(BFsls_ind))
abline(0,1)


# Get decisions based on BF 10

bf01s_dec <- apply(BF01s, c(1,2), function(x) {
    ifelse(x > 10, 0, ifelse(x < 1/10, 1, NA))
})

bf01s_ind_dec <- apply(BF01s_ind, c(1,2), function(x) {
    ifelse(x > 10, 0, ifelse(x < 1/10, 1, NA))
})

# Get decisions based on BF threshold

bf01s_dec <- apply(BFsls, c(1,2), function(x) {
    ifelse(x > 10, 0, ifelse(x < 1/10, 1, NA))
})

bf01s_ind_dec <- apply(BFsls_ind, c(1,2), function(x) {
    ifelse(x > 10, 0, ifelse(x < 1/10, 1, NA))
})

# Decision rate
apply(bf01s_dec, 2, function(x){mean(!is.na(x))})
apply(bf01s_ind_dec, 2, function(x){mean(!is.na(x))})

# Assess decision-accuracy of these things.

truth <- sapply(datasets, function(x){x$params$random_sigma})
truth <- apply(truth, 1:2, function(x) {as.numeric(x > 0)})

apply(truth == bf01s_dec, 2, function(x) {mean(x, na.rm = TRUE)})
apply(truth == bf01s_ind_dec, 2, function(x) {mean(x, na.rm = TRUE)})

