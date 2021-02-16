
test_that("Summary returns desired structures", {
    data(sim_none)
    ds <- sim_none
    fit <- suppressWarnings(mires(myFactor ~ x_1 + x_2 + x_3 + x_4 + x_5 + x_6 + x_7 + x_8, group, ds, iter = 100, cores = 1, chains  = 1))

    sum_fit <- summary(fit)

    expect_s3_class(sum_fit, "summary.mires")
    expect_type(sum_fit$meta, "list")
    expect_type(sum_fit$summary, "list")
    expect_named(sum_fit$summary, c("lambda", "resid", "nu", "resd", "hm", "hm_param", "hm_item", "hm_lambda"))
    for(i in 1:length(sum_fit$summary$lambda[,"Mean"])) {
        expect_gte(sum_fit$summary$lambda$Rhat[i], .5)
    }
})
