
# output <- argument %IfNull% list
`%IfNull%` <- function(object, expr) {
    if(is.null(object)) {
        return(eval(expr))
    } else {
        return(object)
    }
}

nlist <- function(...) {
    mc <- match.call()
    out <- list(...)

    not_named <- is.null(names(out))
    is_named <- if(not_named) {
                    FALSE
                } else {
                    nzchar(names(out))
                }

    args <- as.character(mc)[-1] # Not the fn name.

    if(not_named) {
        names(out) <- args
    } else {
        names(out)[!is_named] <- args[!is_named]
    }

    return(out)

}

prob_to_probs <- function(prob) {
    lower <- (1 - prob) / 2
    upper <- 1 - lower
    probs <- c(lower, upper)
    return(probs)
}

tidy_stanpars <- function(stannames, labs = NULL, ...) {
    dots <- list(...)
    n_replace <- length(dots)
    split <- split_stannames(stannames, labs)
    out <- data.frame(param = split$param,
                      split$indices
                      )
    if(n_replace > 0) {
        cols_to_replace <- names(dots)
        for(n in cols_to_replace) {
            # Replace integers with supplied (ordered) labels in named ...
            out[[n]] <- dots[[n]][out[[n]]]
        }
    }

    return(out)
}

split_stannames <- function(stannames, labs = NULL) {
    param_rex <- r"{(.*)\[(.*)\]}"
    param_names <- gsub(param_rex, "\\1", stannames)
    indices_chr <- gsub(param_rex, "\\2", stannames)
    indices_split <- strsplit(indices_chr, ",")
    indices_split <- lapply(indices_split, as.numeric)
    indices <- do.call(rbind, indices_split)

    if(!is.null(labs)) {
        colnames(indices) <- labs
    } else {
        index_names <- c("row", "col", paste0("arr_", 1:10))
        indices_needed <- ncol(indices)
        colnames(indices) <- index_names[1:indices_needed]
    }

    out <- list(param = param_names, indices = indices)
}

reorder_columns <- function(df, order) {
    cn <- colnames(df)
    col <- 1:ncol(df)
    col_order <- match(order, cn)
    col_removed <- col[-col_order]
    col_new <- c(col_order, col_removed)
    df[, col_new]
}
