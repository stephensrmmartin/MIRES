print.mires <- function(x, ...) {
    .sep()
    .writeLine("MIRES model object")
    .sep()

    .writeLine("Specification:")
    for(f in 1:x$meta$F) {
        print(x$meta$formula[[f]])
    }

    .newline()

    .sep()
    .writeLine("Configuration:")
    .writeLine("\t", "Factors:", ifelse(x$meta$multi, "Multivariate", "Univariate"))
    .writeLine("\t", "Latent Mean and Variance Identification:", ifelse(x$meta$sum_coding, "Sum-to-zero, product-to-one.", "Zero-centered random effects."))
    .writeLine("\t", "Hierarchical inclusion model:", ifelse(x$meta$hmre, "Yes", "No"))
    .writeLine("\t", "Latent Scores Saved:", ifelse(x$meta$save_scores, "Yes", "No"))
}

summary.mires <- function(object, prob = .95, ...) {
    
}

print.summary.mires <- function(x, ...) {
    
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
