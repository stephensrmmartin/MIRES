library(MIRES)

set.seed(14)
J <- 8
K <- 5
n <- 50

fixed <- list(lambda = rep(.7, J), resid_log = rep(log(sqrt(1 - .7^2)), J), nu = rep(0, J))
etadist <- "std"

mipatterns <- list(
    none = list("none"),
    constant = list("constant", .4),
    items = list("items", .4, J/2),
    loadings = list("params", .4, 0),
    resid = list("params", .4, 1),
    nu = list("params", .4, 2)
)

datasets <- lapply(mipatterns, MIRES:::datagen_uni, J = J, K = K, n = n, etadist = etadist, fixed = fixed)
