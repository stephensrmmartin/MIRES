##' @title Print function for mires objects.
##' @param x mires object.
##' @param ... Not used.
##' @return x (Invisibly)
##' @author Stephen R. Martin
##' @export
print.mires <- function(x, ...) {
    .sep()
    .writeLine("MIRES model object")
    .sep()

    .writeLine("Specification:")
    for(f in 1:x$meta$F) {
        .writeLine(.printFormula(x$meta$formula[[f]]))
    }

    .newline()

    .sep()
    .writeLine("Configuration:")
    .writeLine("\t", "Factors:", ifelse(x$meta$multi, "Multivariate", "Univariate"))
    .writeLine("\t", "Latent Mean and Variance Identification:", ifelse(x$meta$sum_coding, "Sum-to-zero, product-to-one.", "Zero-centered random effects."))
    .writeLine("\t", "Hierarchical inclusion model:", ifelse(x$meta$hmre, "Yes", "No"))
    .writeLine("\t", "Latent Scores Saved:", ifelse(x$meta$save_scores, "Yes", "No"))
    if(x$meta$hmre) {
        .writeLine("\t", "Inclusion Model Prior Params:", x$meta$hmre_mu, ", ", x$meta$hmre_scale)
    }

    invisible(x)
}

##' Computes summaries for MIRES objects.
##'
##' Computes summary tables for fixed measurement parameters (loadings, residual SDs, and intercepts) and random effect standard deviations (resd).
##' The printed output includes the posterior mean, median, SD, and .95 (Default) highest density intervals.
##' HDIs were chosen instead of quantile intervals because the random effect SDs can be on the boundary of zero if invariance is plausible.
##' Additionally, other columns exist to help aid decisions about invariance:
##' \describe{
##'   \item{BF01}{Bayes factor of invariance (Variance = 0) to non-invariance (Variance > 0)}
##'   \item{BF10}{Bayes factor of non-invariance (Variance > 0) to invariance (Variance = 0). The inverse of BF01 for convenience}
##'   \item{Pr(SD <= \code{less_than})}{The posterior probability that the random effect SD is less than \code{less_than} (Default: .1). Set \code{less_than} to a value below which you would consider the variance to be effectively ignorable.}
##'   \item{BF(SD <= \code{less_than})}{The Bayes Factor comparing effectively-invariant (SD < \code{less_than}) to non-invariant (SD > \code{less_than}). Set \code{less_than} to a value below which you would consider variance to be effectively ignorable. This uses the encompassing prior approach.}
##' }
##' 
##' @title Summary method for mires object.
##' @param object mires object.
##' @param prob Numeric (Default = .95). Probability mass to be contained in the highest posterior density interval.
##' @param less_than Numeric (Default: .1; positive). Value at which to assess Pr(SD < less_than|D).
##' @param ... Not used.
##' @return summary.mires object. List of meta data and summary. Summary is list of summary tables for all fixed effects parameters.
##' @author Stephen R. Martin
##' @method summary mires
##' @export
summary.mires <- function(object, prob = .95, less_than = .1, ...) {
    dots <- list(...)
    meta <- object$meta
    meta$prob <- prob
    meta$digits <- dots$digits %IfNull% 3

    #################
    # Fixed effects #
    #################

    ## Loadings (Assumes Univariate due to loadings only being a vector, not a matrix!!!)
    lambda <- .summary_table(object,
                            pars = "lambda",
                            prob,
                            labs = "Item",
                            Item = object$meta$indicators)

    ## Resid
    resid <- .summary_table(object,
                            pars = "resid_log",
                            prob,
                            transform = function(x){exp(x)},
                            labs = "Item",
                            Item = object$meta$indicators)
    ## Nu
    nu <- .summary_table(object,
                         pars = "nu",
                         prob,
                         labs = "Item",
                         Item = object$meta$indicators)

    ##########
    # RE-SDs #
    ##########

    # RE-SDs (Assumes Univariate due to loadings only being a vector, not a matrix!!!)
    ## Get loading-codes in row-major order for future proofing: "j_f"
    lambda_enum <- with(object$meta$ind_spec, {
        r <- lapply(1:object$meta$F, function(x){
            paste0(F_ind[x, 1:J_f[x]], "__SEP__", x)
        })
        unlist(r)
    })
    resd_labs <- c(paste0("Loading__SEP__", lambda_enum),
                   paste0("Residual SD__SEP__", 1:object$meta$J),
                   paste0("Intercept__SEP__", 1:object$meta$J)
                   )
    resd <- .summary_table(object,
                           pars = "random_sigma",
                           prob,
                           add_zero = TRUE,
                           labs = "Parameter",
                           Parameter = resd_labs)

    # Munge the __SEP__ values
    param_split <- strsplit(resd[, "Parameter"], "__SEP__")
    resd[, "Parameter"] <- sapply(param_split, function(x) {x[1]})
    resd[, "Item"] <- as.numeric(sapply(param_split, function(x) {x[2]}))
    resd[, "Factor"] <- as.numeric(sapply(param_split, function(x) {x[3]}))
    ## Relabel
    resd[resd$Parameter == "resid_log", "Parameter"] <- "resid"
    resd[, "Item"] <- object$meta$indicators[resd[,"Item"]]
    resd[, "Factor"] <- object$meta$factors[resd[,"Factor"]]
    ## Reorder
    resd <- reorder_columns(resd, c("param", "Parameter", "Item", "Factor"))

    ## Bayes factors
    resd_dfuns <- posterior_density_funs_sigmas(object, add_zero = TRUE)
    bf01 <- sapply(resd_dfuns, function(f) {
        savage_dickey(f, dhmre, x = 0, mu = object$meta$hmre_mu, sigma = object$meta$hmre_scale)
    })
    resd[,"BF01"] <- bf01
    resd[,"BF10"] <- 1 / bf01

    ## Less_than tests
    pr_less_than <- apply(as.matrix(object$fit, pars = "random_sigma"), 2, prob_less_than,
                          less_than = less_than)
    bflts <- apply(as.matrix(object$fit, pars = "random_sigma"), 2, bflt,
                   less_than = less_than,
                   prior_cumul_fun = phmre,
                   mu = object$meta$hmre_mu,
                   sigma = object$meta$hmre_scale
                   )
    resd[,paste0("Pr(SD <= ", less_than, "| D)")] <- pr_less_than
    resd[,paste0("BF(SD <= ", less_than, ")")] <- bflts

    ###############
    # Regularizer #
    ###############

    if(object$meta$hmre) {
        hm <- suppressWarnings(.summary_table(object,
                            pars = "hm_tau",
                            prob,
                            labs = "Parameter")
                            )
        hm_param <- .summary_table(object,
                                pars = "hm_param",
                                prob,
                                labs = "Parameter",
                                Parameter = c("Loading", "Resid", "Intercept"))
        hm_item <- .summary_table(object,
                                pars = "hm_item",
                                prob,
                                labs = "Item",
                                Item = object$meta$indicator)
        hm_lambda <- .summary_table(object,
                                    pars = "hm_lambda",
                                    prob,
                                    labs = "RE"
                                    )
        hm_lambda[,c("Parameter", "Item", "Factor")] <- resd[,c("Parameter", "Item", "Factor")]
        hm_lambda <- reorder_columns(hm_lambda, c("param", "RE", "Parameter", "Item", "Factor"))
    } else {
        hm <- hm_param <- hm_item <- hm_lambda <- NA
    }

    ##########
    # Return #
    ##########

    out <- nlist(meta)
    out$summary <- nlist(lambda,
                 resid,
                 nu,
                 resd,
                 hm,
                 hm_param,
                 hm_item,
                 hm_lambda)

    class(out) <- "summary.mires"
    out
}

##' @title Print method for MIRES summary objects.
##' @param x summary.mires object.
##' @param ... Not used.
##' @return x (Invisibly)
##' @author Stephen Martin
##' @export
print.summary.mires <- function(x, ...) {
    dots <- list(...)
    digits <- dots$digits %IfNull% x$meta$digits

    # Print what print.mires does; only needs meta.
    print.mires(x)
    .sep()

    # Fixed effects
    .writeLine("Fixed Effects")
    .sep()
    ## Loadings (Assumes unidimensional!)
    .writeLine("Loadings")
    .print_sumtab(x$summary$lambda, digits, "param")
    .newline()

    ## Residual SDs
    .writeLine("Residual Standard Deviations (Unlogged Scale)")
    .print_sumtab(x$summary$resid, digits, "param")
    .newline()

    ## Intercepts
    .writeLine("Intercepts")
    .print_sumtab(x$summary$nu, digits, "param")
    .newline()

    ## Random Effect SDs
    .sep()
    .writeLine("Measurement Invariance Assessment")
    .sep()
    .writeLine("Random Effect Standard Deviations (Unlogged Scale)")

    .print_sumtab(x$summary$resd, digits, "param")

    invisible(x)
}


# Internal formatting utilities
.newline <- function(n = 1) {
    cat(rep("\n", n))
}

.sep <- function(sep = "-----") {
    cat(sep)
    .newline()
}

.writeLine <- function(...) {
    dots <- list(...)
    dots <- lapply(dots, as.character)
    do.call(cat, dots)
    ## cat(as.character(x))
    .newline()
}

.printFormula <- function(f) {
    deparse1(f)
}

# Summary-specific helpers

##' @title Compute Highest Posterior Density intervals.
##' @param samps MCMC sample matrix.
##' @param prob Numeric in (0,1).
##' @param add_zero Logical (Default: FALSE) - Whether to add zero to samples. Useful for unidirectional effects.
##' @return Matrix of HDIs.
##' @author Stephen Martin
##' @keywords internal
##' @import HDInterval
.hdi <- function(samps, prob, add_zero = FALSE) {
    credMass <- prob
    if(add_zero) {
        samps <- rbind(0, samps)
    }

    hdis <- t(hdi(samps, credMass = credMass))
    colnames(hdis) <- paste0(c("L", "U"), prob*100)

    return(hdis)
}

.summary_table <- function(mires, pars, prob, add_zero = FALSE, transform = NULL, ...) {
    mcmc <- as.matrix(mires$fit, pars = pars)
    mcmc_array <- as.array(mires$fit, pars = pars)

    if(!is.null(transform)) {
        mcmc <- transform(mcmc)
        mcmc_array <- transform(mcmc_array)
    }

    if(add_zero) {
        mcmc <- rbind(mcmc, 0)
    }

    Mean <- colMeans(mcmc)
    Median <- apply(mcmc, 2, median)
    SD <- apply(mcmc, 2, sd)
    Rhat <- apply(mcmc_array, 3, rstan::Rhat)
    hdis <- .hdi(mcmc, prob, add_zero = add_zero)

    out <- data.frame(Mean, Median, SD, hdis, Rhat)
    out <- cbind(tidy_stanpars(rownames(out), ...), out)
    rownames(out) <- NULL
    return(out)
}

.print_sumtab <- function(x, digits, drop_cols = NULL) {
    which_dc <- match(drop_cols, colnames(x))
    cols <- 1:ncol(x)
    if(length(which_dc) > 0) {
        cols <- (1:ncol(x))[-which_dc]
    }
    x <- x[, cols]

    # Round, because format and print digits are confusing.
    numeric_cols <- sapply(x, is.numeric)
    x[, numeric_cols] <- round(x[, numeric_cols], digits)

    # Print
    print.data.frame(x, row.names = FALSE, na.print = "")
}
