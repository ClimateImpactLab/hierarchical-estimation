## Load the library
source("logspec/logspec.R", chdir=T)

## Just check that everything runs
test.estimate.vcv <- function() {
    ## Load the data
    df <- read.csv("example/true-binned.csv")
    df$log_gdppc <- log(df$gdppc)

    ## Find the coefficients
    result <- estimate.logspec(df$rate, df[, c('bin1', 'bin2', 'bin4', 'bin5')],
                               df[!duplicated(df$adm1), c('meant', 'log_gdppc')],
                               matrix(T, 4, 2), df$adm1, df$adm2)

    ## Determine the VCV for these coefficient
    se.ols <- estimate.vcv(result$betas, result$gammas, result$sigmas,
                           df$rate, df[, c('bin1', 'bin2', 'bin4', 'bin5')],
                           df[!duplicated(df$adm1), c('meant', 'log_gdppc')],
                           matrix(T, 4, 2), df$adm1, df$adm2, use.ols=T)$se

    se.bayes <- estimate.vcv(result$betas, result$gammas, result$sigmas,
                             df$rate, df[, c('bin1', 'bin2', 'bin4', 'bin5')],
                             df[!duplicated(df$adm1), c('meant', 'log_gdppc')],
                             matrix(T, 4, 2), df$adm1, df$adm2, use.ols=F)$se

    ## Apply weights
    estimate.vcv(result$betas, result$gammas, result$sigmas,
                 df$rate, df[, c('bin1', 'bin2', 'bin4', 'bin5')],
                 df[!duplicated(df$adm1), c('meant', 'log_gdppc')],
                 matrix(T, 4, 2), df$adm1, df$adm2, use.ols=T, weights=rexp(nrow(df)))$se
}
