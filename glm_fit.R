fit_glm <- function(expression, time, model) {
    # Function that use linear regression of log transformed
    # expression values to estimate half life of genes.
    # expression: Vector of expression values (numeric)
    # time: vector of time points (numeric)
    # model: either poisson or quasipoisson
    fit <- glm(values ~ time, family = model)
    res <- sum(residuals(fit)^2)
    half.life <- log(2)/(-(coef(fit)[2]))
    return(list(half.life = half.life, residuals = res, fit = fit))
}
