library(RUnit)

## Load the library
source("interpolate/tableapi.R", chdir=T)

tstar <- 22

test.ta.estimate.logspec <- function() {
    ## Load the data
    df <- read.csv("example/true-binned.csv")
    df$log_gdppc <- log(df$gdppc)

    ## Find the coefficients
    result <- ta.estimate.logspec(df, 'rate', 'adm1', 'adm2',
                                  c('bin1', 'bin1', 'bin2', 'bin2', 'bin4', 'bin4', 'bin5', 'bin5'),
                                  c('meant', 'log_gdppc', 'meant', 'log_gdppc', 'meant', 'log_gdppc', 'meant', 'log_gdppc'))
    print(result)

    beta1 <- .1 * (tstar - 13)^2 + .015
    beta2 <- .1 * (tstar - 17.5)^2
    beta4 <- .1 * (tstar - 26.5)^2
    beta5 <- .1 * (tstar - 31)^2 + .015
    checkEqualsNumeric(result$betas, c(beta1, beta2, beta4, beta5), tolerance=.01)

    gamma1 <- -1/200
    gamma2 <- -1/20
    checkEqualsNumeric(result$gammas, rep(c(gamma1, gamma2), 4), tolerance=.02)
}

test.ta.estimate.vcv <- function() {
    ## Load the data
    df <- read.csv("example/true-binned.csv")
    df$log_gdppc <- log(df$gdppc)

    result <- ta.estimate.logspec(df, 'rate', 'adm1', 'adm2',
                                  c('bin1', 'bin1', 'bin2', 'bin2', 'bin4', 'bin4', 'bin5', 'bin5'),
                                  c('meant', 'log_gdppc', 'meant', 'log_gdppc', 'meant', 'log_gdppc', 'meant', 'log_gdppc'))

    ta.estimate.vcv(result$betas, result$gamma, result$sigma,
                    df, 'rate', 'adm1', 'adm2',
                    c('bin1', 'bin1', 'bin2', 'bin2', 'bin4', 'bin4', 'bin5', 'bin5'),
                    c('meant', 'log_gdppc', 'meant', 'log_gdppc', 'meant', 'log_gdppc', 'meant', 'log_gdppc'))
}
