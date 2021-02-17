
test_that("Model configurations change as expected.", {
    data(sim_loadings)
    ds <- sim_loadings
    expect_s3_class(sim_loadings, "data.frame")
    form <- myFactor ~ x_1 + x_2 + x_3 + x_4 + x_5 + x_6 + x_7 + x_8
    iter <- 10
    cores <- 1
    chains <- 1

    config_matrix <- expand.grid(inclusion_model = c("dependent", "independent"),
                                 identification = c("sum_to_zero", "hierarchical"),
                                 save_scores = c(TRUE, FALSE),
                                 stringsAsFactors = FALSE)
    fits <- lapply(1:nrow(config_matrix), function(x) {
        args <- list(formula = form,
                     group = quote(group),
                     data = ds,
                     iter = iter,
                     cores = cores,
                     chains = chains,
                     control = list(adapt_delta = .8))
        args <- c(args, config_matrix[x,])

        suppressWarnings((do.call(mires, args)))
    })

    for(i in 1:length(fits)) {
        # Class
        expect_s3_class(fits[[i]], "mires")
        # Save scores
        expect_equal(fits[[i]]$meta$save_scores, config_matrix[i, "save_scores"])
        # Identification method
        expect_equal(fits[[i]]$meta$sum_coding, config_matrix[i, "identification"] == "sum_to_zero")
        # Inclusion model independent or dependent
        expect_equal(fits[[i]]$meta$hmre, config_matrix[i, "inclusion_model"] == "dependent")
        # All should be multi = FALSE for now.
        expect_false(fits[[i]]$meta$multi)
        expect_false(fits[[i]]$meta$eta_cor_nonmi)
    }
})
