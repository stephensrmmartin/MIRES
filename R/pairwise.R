##' Compute pairwise differences in group-specific measurement parameters.
##' 
##' For a specified set of parameters, this computes all pairwise differences in the random effects across the posterior. Specifically, this computes the posterior differences of groups' parameters, for all parameters. This is useful for comparing groups' estimates under non-invariance.
##'
##' @title Pairwise comparisons of random parameters.
##' @param mires mires object.
##' @param param Character. One of \code{lambda} (loadings), \code{resid} (residual standard deviation on the log scale), or \code{nu} (intercepts).
##' @param prob Numeric (0-1). Probability mass contained within the highest density interval.
##' @param less_than Numeric (Default: .1; positive). Value at which to assess Pr(|difference| < less_than|D).
##' @param groups Character vector (Optional). If specified, will only compute pairwise differences of the specified groups.
##' @param ... Not used.
##' @return Data frame.
##' @author Stephen Martin
##' @export
pairwise <- function(mires, param = c("lambda", "resid", "nu"), prob = .95, less_than = .1, groups = NULL, ...) {
    param <- match.arg(param)
    group_numeric <- seq_len(mires$meta$group$K)
    if(!is.null(groups)) {
        group_numeric <- mires$meta$group$numeric[match(groups, mires$meta$group$data)]
    }

    param_root <- paste0(param, "_random")

    samps <- as.matrix(mires$fit, pars = param_root)

    # Create data.frame for stan param names
    stanpars <- tidy_stanpars(colnames(samps))
    stanpars$par <- colnames(samps)
    # Include only groups in group_numeric
    stanpars <- stanpars[stanpars$row %in% group_numeric, ]
    # Reconstruct the full stan parameter name
    ## stanpars$par <- apply(stanpars, 1, function(x) {
    ##     paste0(x["param"], "[", paste0(as.character(x[2:length(x)]), collapse = ","), "]")
    ## })

    # Create list of stanpars split by non-group parts.
    stanpars_list <- split(stanpars, stanpars[c(-1, -2, -ncol(stanpars))])

    # For each list, use pars of samps in .pairwise_diff_single
    pairwise_list <- lapply(stanpars_list, function(x) {
        .pairwise_diff_single(samps[, x$par])
    })

    # Summarize each
    pairwise_summaries <- lapply(pairwise_list, .summary_table_pairwise,
                                 prob = prob,
                                 hmre_mu = mires$meta$hmre_mu,
                                 hmre_scale = mires$meta$hmre_scale,
                                 less_than = less_than)

    # Add column to each summary for the parameter
    for(i in seq_along(pairwise_list)) {
        pairwise_summaries[[i]]["Param"] <- param
        pairwise_summaries[[i]]["index"] <- names(stanpars_list)[i]
    }

    # Combine
    out <- do.call(rbind, pairwise_summaries)
    rownames(out) <- NULL

    out <- reorder_columns(out, c("Param", "index"))

    # Rename i, j into original group names
    out[,"i"] <- mires$meta$group$data[match(out[,"i"], mires$meta$group$numeric)]
    out[,"j"] <- mires$meta$group$data[match(out[,"j"], mires$meta$group$numeric)]

    out
    
}

.summary_table_pairwise <- function(mcmc, prob, less_than, hmre_mu, hmre_scale, ...) {
    Mean <- colMeans(mcmc)
    Median <- apply(mcmc, 2, median)
    SD <- apply(mcmc, 2, sd)
    hdis <- .hdi(mcmc, prob, add_zero = FALSE)

    # BF01, BF10
    pw_dfuns <- apply(mcmc, 2, logspline::logspline)
    bf01 <- sapply(pw_dfuns, function(l) {logspline::dlogspline(0, l)})
    bf01 <- bf01 / dhmre_pairwise(0, mu = hmre_mu, sigma = hmre_scale)
    bf10 <- 1 / bf01

    # p(Less than), BF(Less than) TODO: BF(Less than) not implemented; needs phmre_pairwise
    plt <- apply(abs(mcmc), 2, prob_less_than, less_than = less_than)


    out <- data.frame(Mean, Median, SD, hdis)
    out[, c("bf01", "bf10")] <- cbind(bf01, bf10)
    out[,paste0("Pr(|i-j|<", less_than, "|D)")] <- plt

    labels <- do.call(rbind, strsplit(rownames(out), "-"))
    colnames(labels) <- c("i", "j")
    out <- cbind(labels, out)

    out
}

##' Computes p_j - p_i not j, for each j.
##'
##' Assumes that mcmc contains params from ONE parameter across the indices whose differences are of interest (e.g, resid_random[1:K, 1]).
##' @title Outer subtraction for given params across MCMC samples.
##' @param mcmc Numeric matrix. MCMC samples.
##' @return MCMC Matrix of all differences. Posterior samples of all possible differences, minus duplicates.
##' @author Stephen R Martin
##' @keywords internal
.pairwise_diff_single <- function(mcmc) {
    stan_indices <- split_stannames(colnames(mcmc))$indices
    diffs <- t(apply(mcmc, 1, .sample_diff))
    if(dim(diffs)[1] == 1) { # Stupid exception for one-comparison case. R drop strikes again.
        diffs <- t(diffs)
    }
    colnames(diffs) <- .sample_diff_labels(stan_indices[,1])
    diffs
}
##' @title Compute all differences of vector.
##' @param x Numeric vector (e.g., one row of mcmc samples)
##' @return Numeric Vector of lower-triangular outer difference matrix.
##' @author Stephen R Martin
##' @keywords internal
.sample_diff <- function(x) {
    out <- outer(x, x, "-")
    out[lower.tri(out)]
}
##' @title Generate labels for all differences of vector.
##' @param mcmcNames 
##' @return Character Vector of lower-triangular outer difference matrix. I.e., labels for .sample_diff.
##' @author Stephen R Martin
##' @keywords internal
.sample_diff_labels <- function(mcmcNames) {
    out <- outer(mcmcNames, mcmcNames, function(x, y) {paste0(x, "-", y)})
    out <- out[lower.tri(out)]
    out
}
