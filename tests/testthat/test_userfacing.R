test_that("User-facing functions return correct structures.", {
    data(sim_loadings)
    ds <- sim_loadings
    expect_s3_class(sim_loadings, "data.frame")

    form <- myFactor ~ x_1 + x_2 + x_3 + x_4 + x_5 + x_6 + x_7 + x_8

    fit <- suppressWarnings(mires(form, group, ds, iter = 100, cores = 1, chains = 1, control=list(adapt_delta = .8)))

    # MIRES object
    expect_s3_class(fit, "mires")
    expect_named(fit, c("meta", "fit", "stan_data"))
    expect_s4_class(fit$fit, "stanfit")
    expect_equal(length(unique(ds$group)), fit$meta$K)
    expect_equal(ncol(ds) - 1, fit$meta$J)

    # Summary
    sum_fit <- summary(fit)
    expect_s3_class(sum_fit, "summary.mires")
    expect_type(sum_fit$meta, "list")
    expect_type(sum_fit$summary, "list")
    expect_named(sum_fit$summary, c("lambda", "resid", "nu", "resd", "hm", "hm_param", "hm_item", "hm_lambda"))
    for(i in 1:length(sum_fit$summary$resd[,"Mean"])) {
        expect_gte(sum_fit$summary$resd$Rhat[i], .5)
    }

    # Ranef
    ranefs <- ranef(fit)
    expect_type(ranefs, "list")
    expect_named(ranefs, c("lambda", "resid", "nu", "eta_mean", "eta_sd"))
    expect_equal(nrow(ranefs$lambda), fit$meta$K * fit$meta$J)

    # Pairwise
    pairwise_out <- pairwise(fit, "nu")
    expect_s3_class(pairwise_out, "data.frame")
    expect_equal(nrow(pairwise_out), (fit$meta$K * (fit$meta$K - 1) / 2) * fit$meta$J)

    pairwise_reduced <- pairwise(fit, param = "lambda", groups = c("1", "2", "3"))
    expect_equal(nrow(pairwise_reduced), (3 * 2 / 2) * fit$meta$J)

})
