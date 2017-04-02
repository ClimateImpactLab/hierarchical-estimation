##setwd("~/research/gcp/hierarchical-estimation/interpolate")

## Load the library
source("logspec.R")

test.estimate.logspec <- function() {
    ## Load the data
    df <- read.csv("../example/true-binned.csv")
    df$log_gdppc <- log(df$gdppc)

    ## Find the coefficients
    result <- estimate.logspec(df$rate, df[, c('bin1', 'bin2', 'bin4', 'bin5')],
                               df[!duplicated(df$adm1), c('meant', 'log_gdppc', 'meant', 'log_gdppc',
                                                          'meant', 'log_gdppc', 'meant', 'log_gdppc')],
                               df$adm1, df$adm2)

    beta1 <- .1 * (tstar - 13)^2
    beta2 <- .1 * (tstar - 17.5)^2
    beta4 <- .1 * (tstar - 26.5)^2
    beta5 <- .1 * (tstar - 31)^2
    checkEquals(result$betas, c(beta1, beta2, beta4, beta5))

    gamma1 <- -1/200
    gamma2 <- -1/20
    checkEquals(result$betas, t(matrix(c(gamma1, gamma2), 2, 4)))
}

## ## Apply weights
## estimate.logspec(df$rate, df[, c('bin1', 'bin2', 'bin4', 'bin5')],
##                  df[!duplicated(df$adm1), c('meant', 'log_gdppc', 'meant', 'log_gdppc',
##                                             'meant', 'log_gdppc', 'meant', 'log_gdppc')],
##                  df$adm1, df$adm2, rexp(nrow(df)))


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

## ## Solve using R optimization method
## estimate.logspec.optim(df$rate, df[, c('bin1', 'bin2', 'bin4', 'bin5')],
##                        df[!duplicated(df$adm2), c('meant', 'log_gdppc', 'meant', 'log_gdppc',
##                                                   'meant', 'log_gdppc', 'meant', 'log_gdppc')],
##                        df$adm2, df$adm2)

## estimate.logspec.optim(df$rate, df[, c('bin1', 'bin2', 'bin4', 'bin5')],
##                        df[!duplicated(df$adm2), c('meant', 'log_gdppc', 'meant', 'log_gdppc',
##                                                   'meant', 'log_gdppc', 'meant', 'log_gdppc')],
##                        df$adm2, df$adm2, initgammas=result$gammas)

## ## Apply weights
## estimate.logspec.optim(df$rate, df[, c('bin1', 'bin2', 'bin4', 'bin5')],
##                        df[!duplicated(df$adm2), c('meant', 'log_gdppc', 'meant', 'log_gdppc',
##                                                   'meant', 'log_gdppc', 'meant', 'log_gdppc')],
##                        df$adm2, df$adm2, rexp(nrow(df)))
