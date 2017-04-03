## Load the library
source("interpolate/logspec.R", chdir=T)

tstar <- 22

test.estimate.logspec <- function() {
    ## Load the data
    df <- read.csv("example/true-binned.csv")
    df$log_gdppc <- log(df$gdppc)

    ## Find the coefficients
    result <- estimate.logspec(df$rate, df[, c('bin1', 'bin2', 'bin4', 'bin5')],
                               df[!duplicated(df$adm1), c('meant', 'log_gdppc', 'meant', 'log_gdppc',
                                                          'meant', 'log_gdppc', 'meant', 'log_gdppc')],
                               df$adm1, df$adm2)

    beta1 <- .1 * (tstar - 13)^2 + .015
    beta2 <- .1 * (tstar - 17.5)^2
    beta4 <- .1 * (tstar - 26.5)^2
    beta5 <- .1 * (tstar - 31)^2 + .015
    checkEqualsNumeric(result$betas, c(beta1, beta2, beta4, beta5), tolerance=.01)

    gamma1 <- -1/200
    gamma2 <- -1/20
    checkEqualsNumeric(result$gammas, t(matrix(c(gamma1, gamma2), 2, 4)), tolerance=.01)

    ## Apply weights
    result2 <- estimate.logspec(df$rate, df[, c('bin1', 'bin2', 'bin4', 'bin5')],
                                df[!duplicated(df$adm1), c('meant', 'log_gdppc', 'meant', 'log_gdppc',
                                                           'meant', 'log_gdppc', 'meant', 'log_gdppc')],
                                df$adm1, df$adm2, rexp(nrow(df)))

    checkEqualsNumeric(result2$betas, c(beta1, beta2, beta4, beta5), tolerance=.02)
    checkEqualsNumeric(result2$gammas, t(matrix(c(gamma1, gamma2), 2, 4)), tolerance=.02)

    ## ## Try to solve without hierarchy
    ## estimate.logspec(df$rate, df[, c('bin1', 'bin2', 'bin4', 'bin5')],
    ##                  df[!duplicated(df$adm2), c('meant', 'log_gdppc', 'meant', 'log_gdppc',
    ##                                             'meant', 'log_gdppc', 'meant', 'log_gdppc')],
    ##                  df$adm2, df$adm2)

    ## ## A cautionary tale: I get the right answer if I provide starting point
    ## estimate.logspec(df$rate, df[, c('bin1', 'bin2', 'bin4', 'bin5')],
    ##                  df[!duplicated(df$adm2), c('meant', 'log_gdppc', 'meant', 'log_gdppc',
    ##                                             'meant', 'log_gdppc', 'meant', 'log_gdppc')],
    ##                  df$adm2, df$adm2, initgammas=result$gammas)
}

test.estimate.logspec.optim <- function() {
    ## Load the data
    df <- read.csv("example/true-binned.csv")
    df$log_gdppc <- log(df$gdppc)

    ## Solve using R optimization method
    result3 <- estimate.logspec.optim(df$rate, df[, c('bin1', 'bin2', 'bin4', 'bin5')],
                                      df[!duplicated(df$adm1), c('meant', 'log_gdppc', 'meant', 'log_gdppc',
                                                                 'meant', 'log_gdppc', 'meant', 'log_gdppc')],
                                      df$adm1, df$adm2)

    beta1 <- .1 * (tstar - 13)^2
    beta2 <- .1 * (tstar - 17.5)^2
    beta4 <- .1 * (tstar - 26.5)^2
    beta5 <- .1 * (tstar - 31)^2
    checkEqualsNumeric(result3$betas, c(beta1, beta2, beta4, beta5), tolerance=1)

    gamma1 <- -1/200
    gamma2 <- -1/20
    checkEqualsNumeric(result3$gammas, t(matrix(c(gamma1, gamma2), 2, 4)), tolerance=1)

    result <- estimate.logspec(df$rate, df[, c('bin1', 'bin2', 'bin4', 'bin5')],
                               df[!duplicated(df$adm1), c('meant', 'log_gdppc', 'meant', 'log_gdppc',
                                                          'meant', 'log_gdppc', 'meant', 'log_gdppc')],
                               df$adm1, df$adm2)

    result4 <- estimate.logspec.optim(df$rate, df[, c('bin1', 'bin2', 'bin4', 'bin5')],
                                      df[!duplicated(df$adm1), c('meant', 'log_gdppc', 'meant', 'log_gdppc',
                                                                 'meant', 'log_gdppc', 'meant', 'log_gdppc')],
                                      df$adm1, df$adm2, initgammas=result$gammas)

    checkEqualsNumeric(result4$betas, c(beta1, beta2, beta4, beta5), tolerance=.01)
    checkEqualsNumeric(result4$gammas, t(matrix(c(gamma1, gamma2), 2, 4)), tolerance=.02)

    ## Apply weights
    result5 <- estimate.logspec.optim(df$rate, df[, c('bin1', 'bin2', 'bin4', 'bin5')],
                                      df[!duplicated(df$adm1), c('meant', 'log_gdppc', 'meant', 'log_gdppc',
                                                                 'meant', 'log_gdppc', 'meant', 'log_gdppc')],
                                      df$adm1, df$adm2, rexp(nrow(df)), initgammas=result$gammas)

    checkEqualsNumeric(result5$betas, c(beta1, beta2, beta4, beta5), tolerance=.02)
    checkEqualsNumeric(result5$gammas, t(matrix(c(gamma1, gamma2), 2, 4)), tolerance=.02)
}
