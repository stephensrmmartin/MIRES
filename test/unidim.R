
library(MIRES)

J <- 10
K <- 20
n <- 40

fixed <- list(lambda = rep(.7, J), resid_log = rep(log(sqrt(1 - .7^2)), J), nu = rep(0, J))
mipattern <- list("items", .4, floor(J/2))
etadist <- "std"

d <- datagen_uni(J, K, n, fixed, mipattern, etadist)
