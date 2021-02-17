
test_that("Random effect outputs are structured correctly.", {
    data(sim_none)
    ds <- sim_none
    fit <- suppressWarnings(mires(myFactor ~ x_1 + x_2 + x_3 + x_4 + x_5 + x_6 + x_7 + x_8, group, ds, cores = 1, chains = 1, iter = 10))

    ranefs <- ranef(fit)
    expect_type(ranefs, "list")
    expect_named(ranefs, c("lambda", "resid", "nu", "eta_mean", "eta_sd"))
    expect_equal(nrow(ranefs$lambda), length(unique(ds$group)) * 8)
})
