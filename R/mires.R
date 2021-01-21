
##' Fits mixed effects measurement models for measurement invariance assessment.
##'
##' MIRES stands for Measurement Invariance assessment with Random Effects and Shrinkage.
##' Specifically, it fits a random-effects measurement model, with the addition of a "Hierarchical Inclusion Model" (the "Shrinkage" part) on the random effect variances.
##' All measurement model parameters (loadings, residual SDs, and intercepts) are modeled as possibly randomly varying across groups.
##' The random effect variances are then hierarchically modeled from a half-normal distribution with location zero, and estimated scale.
##' The prior scale controls whether the RE variance is effectively zero (invariant) or not (non-invariant).
##' The scale is itself hierarchically modeled and partially pooled.
##' Therefore, whether a parameter is invariant (the variance being effectively zero) or not (the variance permitted to be non-zero) is informed by all other parameters.
##' Currently, we assume that each parameter informs the invariance of other similar parameters (presence of variance in loadings informs the presence of variance in other loadings), and of similar items (non-invariance of item j parameters informs other parameters for item j).
##' The benefit of this is that information about the presence or absence of invariance of parameters is increased, allowing for more certain decisions about the presence and magnitude of invariance.
##' This is in contrast to the typical random effect approach, wherein it is assumed that all RE variances are independent of one another.
##' @title Fit mixed effects measurement model for invariance assessment.
##' @param formula Formula (or list of formulas). LHS is the factor name, and RHS contains indicators.
##' @param group Grouping variable (as raw name). Grouping variable over which to assess invariance.
##' @param data data.frame. Must contain the indicators specified in formula, and the grouping variable.
##' @param ... Further arguments to \code{\link[rstan]{sampling}}.
##' @return mires object.
##' @import rstan
##' @importFrom parallel detectCores
##' @author Stephen R. Martin
##' @export
mires <- function(formula, group, data, ...) {
    dots <- list(...)

    # Get parsed data structures
    d <- .parse_formula(formula, group, data)

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
    multi <- d$meta$F > 1
    sum_coding <- dots$sum_coding %IfNull% TRUE
    eta_cor_nonmi <- dots$eta_cor_nonmi %IfNull% FALSE
    prior_only <- dots$prior_only %IfNull% FALSE
    save_scores <- dots$save_scores %IfNull% FALSE
    hmre <- dots$hmre %IfNull% TRUE
    hmre_mu <- dots$hmre_mu %IfNull% 0.0
    hmre_scale <- dots$hmre_scale %IfNull% .25
    dots[c("sum_coding",
        "eta_cor_nonmi",
        "prior_only",
        "save_scores",
        "hmre",
        "hmre_mu",
        "hmre_scale")] <- NULL

    # Save config options to metadata list
    d$meta <- c(d$meta, nlist(
                            multi,
                            sum_coding,
                            eta_cor_nonmi,
                            prior_only,
                            save_scores,
                            hmre,
                            hmre_mu,
                            hmre_scale
                        )
                )

    stan_args$data$eta_cor_nonmi <- eta_cor_nonmi
    stan_args$data$prior_only <- prior_only
    stan_args$data$hmre_mu <- hmre_mu
    stan_args$data$hmre_scale <- hmre_scale

    ## Select model
    if(multi) { # Multidimensional
        model <- "redifhm_multi_hier"
    } else if(sum_coding) {
        model <- "redifhm_sum"
    } else { # Fallback
        model <- "redifhm_hier" 
    }
    stan_args$object <- stanmodels[[model]] 

    ## Select params
    ### Shared params
    pars <- c("lambda", "resid_log", "nu",
              "lambda_random", "resid_random", "nu_random",
              "eta_mean", "eta_sd",
              "hm_tau", "hm_param", "hm_item", "hm_lambda")
    if(save_scores) {
        pars <- c(pars, "eta")
    }
    pars <- c(pars, "RE_cor", "random_sigma")

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
    group <- substitute(group)
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
