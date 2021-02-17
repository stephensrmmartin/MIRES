
test_that("Group argument is parsed correctly.", {
    data <- data.frame(x1 = rnorm(10), x2 = rnorm(10), aGroupName = rep(1:2, each = 5))
    form <- fac ~ x1 + x2
    f <- function(f, group, data) {
        .parse_formula(f, substitute(group), data)
    }
    parsed <- f(form, aGroupName, data)
    ## parsed <- .parse_formula(form, aGroupName, data)
    expect_equal(parsed$meta$group$name, "aGroupName")
    expect_s3_class(parsed$meta$group$data, "factor")
})






