
##' Fits mixed effects measurement models for measurement invariance assessment.
##'
##' MIRES stands for Measurement Invariance assessment with Random Effects and Shrinkage.
##' Unlike other measurement invariance approaches, the MIRES model assumes all measurement model parameters (loadings, residual SDs, and intercepts) can randomly vary across groups --- It is a mixed effects model on all parameters.
##' Unlike most mixed effects models, the random effect variances are themselves also hierarchically modeled from a half-normal distribution with location zero, and a scaling parameter.
##' This scaling parameter allows for rapid shrinkage of variance toward zero (invariance), while allowing variance if deemed necessary (non-invariance).
##' 
##' The scaling parameter (an estimated quantity) controls whether the RE variance is effectively zero (invariant) or not (non-invariant).
##' Therefore, the random effect variances are regularized.
##' When \code{inclusion_model} is \code{dependent} (Default), the scaling parameters are hierarchically modeled.
##' By doing so, the invariance or non-invariance of a parameter is informed by other parameters with shared characteristics.
##' Currently, we assume that each parameter informs the invariance of other similar parameters (presence of variance in some loadings increases the probability of variance in other loadings), and of similar items (non-invariance of item j parameters informs other parameters for item j).
##' This approach increases the information available about presence or absence of invariance, allowing for more certain decisions.
##'
##' The "Hierarchical inclusion model" on the random effect variance manifests as a hierarchical prior.
##' When a dependent inclusion model is specified, then the hierarchical prior on random effect SDs is:
##' \deqn{p(\sigma_p | \exp(\tau)) = \mathcal{N}^+(\sigma_p | 0, \exp(\tau))}
##' \deqn{\tau = \tau_c + \tau_{param} + \tau_{item} + \tau_p}
##' \deqn{\tau_* \sim \mathcal{N}(\mu_h, \sigma_h)}
##' Therefore, the regularization of each RE-SD is shared between all RE-SDs (tau_c), all RE-SDs of the same parameter type (tau_param), and all RE-SDs of the same item (tau_item).
##' 
##' When an independent inclusion model is specified (\code{inclusion_model} is "independent"), only the independent regularization term \eqn{\tau_p} is included.
##' The prior is then scaled so that the marginal prior on each \eqn{p(\sigma_p)} remains the same.
##' In this case, RE-SDs cannot share regularization intensities between one another.
##'
##' The inclusion model hyper parameters (mu_h, sigma_h) can be specified, but we recommend the default as a relatively sane, but weakly informative prior.
##' @title Fit mixed effects measurement model for invariance assessment.
##' @param formula Formula. LHS is the factor name, and RHS contains indicators.
##' @param group Grouping variable (symbol). Grouping variable over which to assess invariance.
##' @param data data.frame. Must contain the indicators specified in formula, and the grouping variable.
##' @param inclusion_model String (Default: dependent). If dependent, then the regularization of RE-SDs are dependent (See Details). If independent, then regularization is per-parameter. This is useful for comparing a dependent inclusion model to a non-dependent inclusion model. Note that adaptive regularization occurs regardless (until a non-regularized version is implemented).
##' @param identification String (Default: sum_to_zero). If \code{hierarchical}, then latent means and (log) SDs are identified as zero-centered random effects. If \code{sum_to_zero}, then latent means are identified by a sum-to-zero constraint, and (log) latent SDs are identified by a sum-to-zero constraint.
##' @param save_scores Logical (Default: FALSE). If TRUE, latent scores for each observation are estimated. If FALSE (Default), latent scores are marginalized out; this can result in more efficient sampling and faster fits, due to the drastic reduction in estimated parameters. Note that the random effects for each group are always estimated, and are not marginalized out.
##' @param prior_only Logical (Default: FALSE). If TRUE, samples are drawn from the prior.
##' @param prior Numeric vector (Default: c(0, .25)). The location and scale parameters for the hierarchical inclusion model.
##' @param ... Further arguments to \code{\link[rstan]{sampling}}.
##' @return mires object.
##' @import rstan
##' @importFrom parallel detectCores
##' @import stats
##' @author Stephen R. Martin
##' @export
##' @references
##' Stan Development Team (2020). RStan: the R interface to Stan. R package version 2.21.1. https://mc-stan.org
##'
##' Martin, S. R., Williams, D. R., & Rast, P. (2019, June 18). Measurement Invariance Assessment with Bayesian Hierarchical Inclusion Modeling. <doi:10.31234/osf.io/qbdjt>
##' @examples
##' \donttest{
##' data(sim_loadings) # Load simulated data set
##' head(sim_loadings) # 8 indicators, grouping variable is called "group"
##' 
##' # Fit MIRES to simulated data example.
##' # Assume factor name is, e.g., agreeableness.
##' fit <- mires(agreeableness ~ x_1 + x_2 + x_3 + x_4 + x_5 + x_6 + x_7 + x_8,
##'              group = group,
##'              data = sim_loadings, chains = 2, cores = 2)
##'
##' # Summarize fit
##' summary(fit)
##'
##' # Compare all groups' loadings:
##' pairwise(fit, param = "lambda")
##' # Compare groups "2" and "3" only:
##' pairwise(fit, param = "lambda", groups = c("2", "3"))
##'
##' # Get random effects:
##' fit_ranefs <- ranef(fit)
##' # Look at random effects of loadings:
##' fit_ranefs$lambda
##' }
mires <- function(formula,
                  group,
                  data,
                  inclusion_model = c("dependent", "independent"),
                  identification = c("sum_to_zero", "hierarchical"),
                  save_scores = FALSE,
                  prior_only = FALSE,
                  prior = c(0, .25),
                  ...) {
    dots <- list(...)

    # Get parsed data structures
    d <- .parse_formula(formula, substitute(group), data)

    # Initialize Stan arguments
    stan_args <- list(data = d$stan_data,
                      cores = detectCores(),
                      chains = dots$chains %IfNull% 4,
                      control = dots$control %IfNull% list(adapt_delta = .95),
                      init = dots$init %IfNull% 0
                      )
    stan_args$control$adapt_delta <- stan_args$control$adapt_delta %IfNull% list(adapt_delta = .95)
    # Avoid arg duplicates
    dots[names(dots) %in% names(stan_args)] <- NULL

    # Model Configuration
    inclusion_model <- match.arg(inclusion_model)
    ident <- match.arg(identification)
    hmre <- inclusion_model == "dependent"
    save_scores <- save_scores
    marginalize <- !save_scores
    prior_only <- prior_only
    hmre_mu <- prior[1]
    hmre_scale <- prior[2]
    sum_coding <- ident == "sum_to_zero"

    ## For multidimensional models (Not implemented yet)
    multi <- d$meta$F > 1
    eta_cor_nonmi <- dots$eta_cor_nonmi %IfNull% FALSE # Allow latent correlations to vary (TRUE) or not (FALSE)
    dots[c(
        "eta_cor_nonmi"
        )] <- NULL # Remove these from dots.

    # Save config options to metadata list
    d$meta <- c(d$meta, nlist(
                            multi,
                            eta_cor_nonmi,
                            sum_coding,
                            prior_only,
                            save_scores,
                            hmre,
                            hmre_mu,
                            hmre_scale,
                            marginalize
                        )
                )

    # Push config options to Stan
    stan_args$data <- within(stan_args$data,{
        eta_cor_nonmi = eta_cor_nonmi
        prior_only = prior_only
        hmre_mu = hmre_mu
        hmre_scale = hmre_scale
        use_hmre = hmre
        marginalize = marginalize
        sum_coding = sum_coding
    })

    ## Select model
    stan_args$object <- stanmodels[["redifhm_all"]]

    ## Select params
    ### Shared params
    pars <- c("lambda", "resid_log", "nu",
              "lambda_random", "resid_random", "nu_random",
              "eta_mean", "eta_sd", "RE_cor", "random_sigma")
    if(hmre) {
        pars <- c(pars, "hm_tau", "hm_param", "hm_item", "hm_lambda")
    }
    if(save_scores) {
        pars <- c(pars, "eta")
    }
    if(!sum_coding) {
        pars <- c(pars, "eta_random_sigma")
    }

    stan_args$pars <- pars

    stanOut <- do.call(sampling, c(stan_args, dots))

    out <- list(meta = d$meta,
                fit = stanOut,
                stan_data = d$stan_data)

    class(out) <- "mires"
    return(out)
}

##' @title Parse formula (list).
##' @param formula Formula or list of formulas.
##' @param group Raw name for group.
##' @param data data.frame.
##' @return List containing meta-data and stan data.
##' @import Formula
##' @author Stephen R. Martin
##' @keywords internal
.parse_formula <- function(formula, group, data) {

    formList <- formula

    # Convert to list if not already
    if(!is.list(formList)) {
        formList <- list(formList)
    }
    formList <- lapply(formList, as.Formula)

    # Get group
    group <- group
    group_string <- deparse(group)

    # Get combined formula for model.frame
    RHS <- .combine_RHS(formList)

    # Get relevant variable names
    factors <- .formula_names(formList)$factor
    vars <- c(group_string, unique(do.call(c, .formula_names(formList, terms = FALSE)$indicator)))

    # Reduce data to relevant vars
    data.sub <- data[,vars]
    ## Alt:
    ## data.sub <- get_all_vars(RHS, data, group = group) ## But creates column named group, with group data in it; vs a column named group_string.

    # Remove missings
    data.complete <- data.sub[complete.cases(data.sub),]
    n.removed <- nrow(data.sub) - nrow(data.complete)
    if(n.removed > 0) {
        warning("Removed", n.removed, "incomplete cases.")
    }

    # Restructure group data
    group_data <- as.factor(data.complete[, group_string])
    group_data_numeric <- as.numeric(group_data)
    group_k <- length(unique(group_data_numeric))

    group <- list(
        name = group_string,
        data = group_data,
        numeric = group_data_numeric,
        K = group_k
    )

    # Model matrix
    mm <- model.matrix(RHS, data.complete)[, -1] # Remove intercept

    # Indicator specification
    ind_spec <- .indicator_spec(formList, mm)

    # Meta-data
    N <- nrow(data.complete)
    J <- ncol(mm)
    F <- length(formList)
    K <- group$K

    meta <- list(N = N,
                 J = J,
                 F = F,
                 K = K,
                 group = group,
                 factors = do.call(c, factors),
                 indicators = colnames(mm),
                 ind_spec = ind_spec,
                 data = data.complete,
                 formula = formList)

    stan_data <- list(N = N,
                      J = J,
                      F = F,
                      K = K,
                      group = group$numeric,
                      x = mm,
                      J_f = ind_spec$J_f,
                      F_ind = ind_spec$F_ind
                      )

    out <- list(meta = meta,
                stan_data = stan_data
                )
    return(out)
}

##' @title Get terms from formula list
##' @param formList 
##' @param terms 
##' @return List containing factor (factor names) and indicator (indicator names) as lists.
##' @author Stephen R. Martin
##' @keywords internal
.formula_names <- function(formList, terms = TRUE) {
    fname <- lapply(formList, .formula_lhs)
    iname <- lapply(formList, .formula_rhs)
    out <- list(factor = fname, indicator = iname)

    return(out)
}
##' @title Get the one-length LHS of formula as string.
##' @param formula Formula.
##' @return String. LHS variable of formula. Formula must have only one name on LHS.
##' @author Stephen R. Martin
##' @keywords internal
.formula_lhs <- function(formula) {
    all.vars(formula)[1]
}

##' @title Get RHS of formula as character vector.
##' @param formula Formula.
##' @param terms Logical. Whether to return the formula expressions, or variables (FALSE). I.e., "I(x^2)" instead of "x".
##' @return Character vector.
##' @author Stephen R. Martin
##' @keywords internal
.formula_rhs <- function(formula, terms = TRUE) {
    if(terms) {
        labels(terms(formula))
    } else {
        all.vars(formula)[-1]
    }
}

##' @title Combine all unique RHS entries into one RHS formula.
##' @param formList List of formulas.
##' @return Formula. RHS only.
##' @author Stephen R. Martin
##' @import Formula
##' @keywords internal
.combine_RHS <- function(formList) {
    formula_names <- .formula_names(formList, terms = TRUE)

    iname <- formula_names$indicator
    iname <- do.call(c, iname)
    iname <- unique(iname)

    rhs <- paste0("~", paste(iname, collapse = "+"))
    rhs <- as.Formula(rhs)

    return(rhs)
}

##' Generates the "indicator spec" used by Stan
##'
##' The indicator spec consists of two parts.
##' The first part is J_f, or the number of indicators under each factor.
##' The second part is an [F, J] array, wherein each row defines the 1:J_f[f] columns of the indicator matrix belonging to the factor.
##' Example:
##' [1, 3, 5, 0, 0, 0]: J_f[1] = 3
##' [2, 4, 6, 0, 0, 0]: J_f[2] = 3
##' [1, 2, 3, 4, 0, 0]: J_f[3] = 4; J = 6; F = 3
##' @title Generates indicator spec list.
##' @param formList 
##' @param mm 
##' @return List containing J_f (Indicators per factor; numeric vector) and F_ind (FxJ Numeric Matrix, where F_ind[f,1:J_f] gives the column indices of the model matrix corresponding to factor f.)
##' @author Stephen R. Martin
##' @keywords internal
.indicator_spec <- function(formList, mm) {
    F <- length(formList)
    J <- ncol(mm)

    J_f <- array(0, F)
    F_ind <- matrix(0, F, J)

    formNames <- .formula_names(formList, terms = TRUE)
    

    for(f in 1:F) {
        J_f[f] <- length(formNames$indicator[[f]])
        F_ind[f, 1:J_f[f]] <- match(formNames$indicator[[f]], colnames(mm))
    }

    out <- list(J_f = J_f, F_ind = F_ind)
    return(out)
    
}
