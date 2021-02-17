
test_that("Pairwise returning correct structures.", {
    data(sim_intercepts)
    ds <- sim_intercepts
    K <- length(unique(ds$group))
    K_unique <- K*(K-1) / 2
    fit <- suppressWarnings(mires(myFactor ~ x_1 + x_2 + x_3 + x_4 + x_5 + x_6 + x_7 + x_8, group, ds, iter = 100, cores = 1, chains  = 1))

    expect_s3_class(fit, "mires")
    pairwise_out <- pairwise(fit, "nu")
    expect_s3_class(pairwise_out, "data.frame")
    expect_equal(nrow(pairwise_out), K_unique * 8)

    pairwise_reduced <- pairwise(fit, param = "lambda", groups = c("1", "2", "3"))
    expect_equal(nrow(pairwise_reduced), 3 * 2 / 2 * 8)
})
