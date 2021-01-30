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
etadist <- "std"

## mipattern <- list("none")

d <- MIRES:::datagen_uni(J, K, n, fixed, mipattern, etadist = etadist)

ds <- d$df

# HMRE
fit <- mires(myfactor ~ x_1 + x_2 + x_3 + x_4 + x_5 + x_6 + x_7 + x_8 + x_9 + x_10, group = group, ds, iter = 1000, prior_only = FALSE, hmre_mu = 0, hmre_scale = .25)

# Hierarchical HMRE
fit_hier <- mires(myfactor ~ x_1 + x_2 + x_3 + x_4 + x_5 + x_6 + x_7 + x_8 + x_9 + x_10, group = group, ds, iter = 1000, prior_only = FALSE, hmre_mu = 0, hmre_scale = .25, sum_coding = FALSE)

# Hierarchical non-HMRE
fit_hier_re <- mires(myfactor ~ x_1 + x_2 + x_3 + x_4 + x_5 + x_6 + x_7 + x_8 + x_9 + x_10, group = group, ds, iter = 1000, prior_only = FALSE, hmre_mu = 0, hmre_scale = .25, sum_coding = FALSE, hmre = FALSE)

# Hierarchical Marg.
fit_hier_marg <- mires(myfactor ~ x_1 + x_2 + x_3 + x_4 + x_5 + x_6 + x_7 + x_8 + x_9 + x_10, group = group, ds, iter = 1000, prior_only = FALSE, hmre_mu = 0, hmre_scale = .25, sum_coding = FALSE, marginalize=TRUE)

##################
# COMBINED MODEL #
##################

# HMRE (Combined)
fit <- mires(myfactor ~ x_1 + x_2 + x_3 + x_4 + x_5 + x_6 + x_7 + x_8 + x_9 + x_10, group = group, ds, iter = 1000, prior_only = FALSE, hmre_mu = 0, hmre_scale = .25, combined = TRUE)

# Hierarchical HMRE (Combined)
fit_hier <- mires(myfactor ~ x_1 + x_2 + x_3 + x_4 + x_5 + x_6 + x_7 + x_8 + x_9 + x_10, group = group, ds, iter = 1000, prior_only = FALSE, hmre_mu = 0, hmre_scale = .25, sum_coding = FALSE, combined = TRUE)

# Hierarchical non-HMRE (Combined)
fit_hier_re <- mires(myfactor ~ x_1 + x_2 + x_3 + x_4 + x_5 + x_6 + x_7 + x_8 + x_9 + x_10, group = group, ds, iter = 1000, prior_only = FALSE, hmre_mu = 0, hmre_scale = .25, sum_coding = FALSE, hmre = FALSE, combined = TRUE)

# Hierarchical Marg. (Combined)
fit_hier_marg <- mires(myfactor ~ x_1 + x_2 + x_3 + x_4 + x_5 + x_6 + x_7 + x_8 + x_9 + x_10, group = group, ds, iter = 1000, prior_only = FALSE, hmre_mu = 0, hmre_scale = .25, sum_coding = FALSE, marginalize=TRUE, combined = TRUE)


########
# Misc #
########

# Get mean times (in minutes)
mean(rowMeans(rstan::get_elapsed_time(fit_hier_marg$fit))) / 60
mean(rowMeans(rstan::get_elapsed_time(fit_hier$fit))) / 60
mean(rowMeans(rstan::get_elapsed_time(fit_hier_re$fit))) / 60


# Plot RESDs diffs
plot(1:30, y = summary(fit_hier_marg)$summary$resd[,"Mean"], ylim = c(0, .5), col = "blue")
points(1:30, y = summary(fit_hier)$summary$resd[,"Mean"], ylim = c(0, .2), col = "red")
points(1:30, y = summary(fit_hier_re)$summary$resd[,"Mean"], ylim = c(0, .2), col = "green")
legend("topright",legend = c("Marg", "Non-Marg", "RE [Non-Marg]"), col = c("blue", "red", "green"), pch = 1)
## CrI
segments(x0 = 1:30,
         y0 = summary(fit_hier_marg)$summary$resd[,"L95"],
         y1 = summary(fit_hier_marg)$summary$resd[,"U95"],
         col = "blue")
segments(x0 = 1:30,
         y0 = summary(fit_hier)$summary$resd[,"L95"],
         y1 = summary(fit_hier)$summary$resd[,"U95"],
         col = "red")
segments(x0 = 1:30,
         y0 = summary(fit_hier_re)$summary$resd[,"L95"],
         y1 = summary(fit_hier_re)$summary$resd[,"U95"],
         col = "green")


### BF01s
plot(1:30, y = summary(fit_hier_marg)$summary$resd[,"BF01"], col = "blue")
points(1:30, y = summary(fit_hier)$summary$resd[,"BF01"], col = "red")
points(1:30, y = summary(fit_hier_re)$summary$resd[,"BF01"], col = "green")
legend("topright",legend = c("Marg", "Non-Marg", "RE [Non-Marg]"), col = c("blue", "red", "green"), pch = 1)
