library(MIRES)

set.seed(14)
J <- 8
K <- 5
n <- 50

fixed <- list(lambda = rep(.7, J), resid_log = rep(log(sqrt(1 - .7^2)), J), nu = rep(0, J))
etadist <- NULL

mipatterns <- list(
    none = list("none"),
    all = list("constant", .4),
    items = list("items", .4, J/2),
    loadings = list("params", .4, 0),
    resid = list("params", .4, 1),
    intercepts = list("params", .4, 2)
)

datasets <- lapply(mipatterns, MIRES:::datagen_uni, J = J, K = K, n = n, etadist = etadist, fixed = fixed)

names(datasets) <- paste0("sim_", names(datasets))

for(d in names(datasets)) {
    assign(d, datasets[[d]]$df)
}

usethis::use_data(sim_none)
usethis::use_data(sim_all)
usethis::use_data(sim_items)
usethis::use_data(sim_intercepts)
usethis::use_data(sim_loadings)
usethis::use_data(sim_resid)
