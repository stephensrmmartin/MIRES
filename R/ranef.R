##' @title Extract random effects of each group from MIRES model.
##' @param object mires object.
##' @param prob Numeric (Default: .95). Amount of probability mass to contain within the credible interval.
##' @param ... Not used.
##' @return List containing summaries of lambda, (log) residual SDs, nu, latent mean, and (log) latent SD random effects.
##' @author Stephen R Martin
##' @importFrom nlme ranef
##' @export ranef
##' @export
##' @method ranef mires
##' @aliases ranef
ranef.mires <- function(object, prob = .95, ...) {
    group_name <- object$meta$group$name

    ## Loadings (Assumes Univariate due to loadings only being a vector, not a matrix!!!)
    lambda_random <- .summary_table(object,
                             pars = "lambda_random",
                             prob,
                             labs = c(group_name, "Item"),
                             Item = object$meta$indicators)

    ## Resid_log
    resid_random <- .summary_table(object,
                             pars = "resid_random",
                             prob,
                             labs = c(group_name, "Item"),
                             Item = object$meta$indicators)

    ## Nu (Intercepts)
    nu_random <- .summary_table(object,
                             pars = "nu_random",
                             prob,
                             labs = c(group_name, "Item"),
                             Item = object$meta$indicators)

    ## Eta mean (Assumes Univariate due to eta_mean only being a vector, not a matrix!!!)
    eta_mean <- .summary_table(object,
                             pars = "eta_mean",
                             prob,
                             labs = group_name
                             )
    eta_mean[, group_name] <- levels(object$meta$group$data)[eta_mean[, group_name]]

    ## Eta SD (Assumes Univariate due to eta_mean only being a vector, not a matrix!!!)
    eta_sd <- .summary_table(object,
                             pars = "eta_sd",
                             prob,
                             labs = group_name)
    eta_sd[, group_name] <- levels(object$meta$group$data)[eta_sd[, group_name]]

    out <- list(lambda = lambda_random,
                resid = resid_random,
                nu = nu_random,
                eta_mean = eta_mean,
                eta_sd = eta_sd
                )
    out
}
