library(RUnit)

## Load the library
source("logspec/logspec.R", chdir=T)

tstar <- 22

test.gradient <- function() {
    df <- read.csv("../example/true-binned.csv")
    df$log_gdppc <- log(df$gdppc)

    result1 <- estimate.logspec(df$rate, df[, c('bin1', 'bin2', 'bin4', 'bin5')],
                                df[!duplicated(df$adm1), c('meant', 'log_gdppc')],
                                matrix(T, 4, 2), df$adm1, df$adm2)

    result2 <- estimate.logspec.optim(df$rate, df[, c('bin1', 'bin2', 'bin4', 'bin5')],
                                      df[!duplicated(df$adm1), c('meant', 'log_gdppc')],
                                      matrix(T, 4, 2), df$adm1, df$adm2)

    gradient1 <- get.gamma.gradient(df$rate, df[, c('bin1', 'bin2', 'bin4', 'bin5')],
                                    df[!duplicated(df$adm1), c('meant', 'log_gdppc')],
                                    matrix(T, 4, 2), df$adm1, df$adm2,
                                    result1$betas, result1$gammas, result1$sigma)

    gradient2 <- get.gamma.gradient(df$rate, df[, c('bin1', 'bin2', 'bin4', 'bin5')],
                                    df[!duplicated(df$adm1), c('meant', 'log_gdppc')],
                                    matrix(T, 4, 2), df$adm1, df$adm2,
                                    result2$betas, result2$gammas, result2$sigma)

    checkTrue(all(abs(gradient1) < abs(gradient2)))
}

test.without.gradient <- function() {
    df <- read.csv("../example/true-binned.csv")
    df$log_gdppc <- log(df$gdppc)

    result1 <- estimate.logspec(df$rate, df[, c('bin1', 'bin2', 'bin4', 'bin5')],
                                df[!duplicated(df$adm1), c('meant', 'log_gdppc')],
                                matrix(T, 4, 2), df$adm1, df$adm2)

    print("Optimizing without gradient.")

    result1 <- estimate.logspec(df$rate, df[, c('bin1', 'bin2', 'bin4', 'bin5')],
                                df[!duplicated(df$adm1), c('meant', 'log_gdppc')],
                                matrix(T, 4, 2), df$adm1, df$adm2)

    ptm1 <- proc.time()

    optim1 <- estimate.logspec.gammaoptim.nograd(df$rate, df[, c('bin1', 'bin2', 'bin4', 'bin5')],
                                          df[!duplicated(df$adm1), c('meant', 'log_gdppc')],
                                          matrix(T, 4, 2), df$adm1, df$adm2, result1$sigmas)

    time1 <- proc.time() - ptm1

    print("Optimizing with gradient.")

    ptm2 <- proc.time()

    optim2 <- estimate.logspec.gammaoptim(df$rate, df[, c('bin1', 'bin2', 'bin4', 'bin5')],
                                           df[!duplicated(df$adm1), c('meant', 'log_gdppc')],
                                           matrix(T, 4, 2), df$adm1, df$adm2, result1$sigmas)

    time2 <- proc.time() - ptm

    checkTrue(all(time2 < time1))
}
