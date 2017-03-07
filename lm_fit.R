fit_lm <- function(expression, time) {
    # Function that use linear regression of log transformed
    # expression values to estimate half life of genes.
    # expression: Vector of expression values (numeric)
    # time: vector of time points (numeric)
    fit <- lm(log(expression) ~ time)
    res <- sum(residuals(fit)^2)
    half.life <- log(2)/(-(coef(fit)[2]))
    return(list(half.life = half.life, residuals = res, fit = fit))
}
