
test_that("posterior_density_funs_sigmas works as intended with sufficient MCMC samples", {
    
          mcmc_samples <- matrix(abs(rnorm(200)), ncol = 2)
          colnames(mcmc_samples) <- paste0("random_sigma", "[", 1:ncol(mcmc_samples),"]")
          test_obj <- list(fit = mcmc_samples)
          out <- posterior_density_funs_sigmas(test_obj)
          for(i in 1:length(out)) {
              expect_type(out[[i]], "closure")
          }
}
)


test_that("posterior_density_funs_sigmas calls DP when needed.", {
    mcmc_samples <- matrix(c(.01, .05, 1), nrow = 90, ncol = 2)
    colnames(mcmc_samples) <- paste0("random_sigma", "[", 1:ncol(mcmc_samples),"]")
    test_obj <- list(fit = mcmc_samples)
    # Not needing a good fit; just a quick one
    out <- suppressWarnings(posterior_density_funs_sigmas(test_obj, tol_rel_obj = 1))
    for(i in 1:length(out)) {
        expect_type(out[[i]], "closure")
    }
})


