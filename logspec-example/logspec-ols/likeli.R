##setwd("~/research/gcp/hierarchical-estimation/logspec-ols")

do.tests <- F

calc.expected <- function(stan.data, alphas, betas, gammas) {
    obsmean <- alphas[stan.data$nn]
    for (kk in 1:stan.data$K)
        obsmean <- obsmean + betas[kk] * stan.data$xx[kk, ] * exp(stan.data$zz[kk, , ] %*% gammas[kk, ])

    as.numeric(obsmean)
}

calc.likeli <- function(stan.data, alphas, betas, gammas, sigma) {
    obsmean <- calc.expected(stan.data, alphas, betas, gammas)
    sum(dnorm(stan.data$yy, obsmean, sigma[stan.data$mm], log=T))
}

calc.likeli.partial <- function(stan.data, betas, gammas, sigma) {
    obsmean <- 0
    for (kk in 1:stan.data$K) {
        dmxxexzz <- regional.demean(stan.data$xx[kk, ] * exp(stan.data$zz[kk, , ] %*% gammas[kk, ]), stan.data$nn)
        obsmean <- obsmean + betas[kk] * dmxxexzz
    }

    sum(dnorm(regional.demean(stan.data$yy, stan.data$nn), obsmean, sigma[stan.data$mm], log=T))
}

if (do.tests) {
    source("../example/logspec-data.R", chdir=T)

    teff0 <- .1 # Size of temperature effect
    tstar <- 21 # Lowest-effect temperature
    betas <- teff0 * (c(12, 16.5, 25.5, 30) - tstar)^2 - teff0 * (21 - tstar)^2
    gammas <- t(matrix(rep(c(-.005, -.05), 4), 2, 4))

    df$exprate <- calc.expected(stan.data, df$const[!duplicated(df$adm2)], betas, gammas) + df$dropeff
    ##df$exprate <- df$const + df$dropeff + betas[1] * exp(-df$log_gdppc / 20 - df$meant / 200) * df$bin1 + betas[2] * exp(-df$log_gdppc / 20 - df$meant / 200) * df$bin2 + betas[3] * exp(-df$log_gdppc / 20 - df$meant / 200) * df$bin4 + betas[4] * exp(-df$log_gdppc / 20 - df$meant / 200) * df$bin5

    print(calc.likeli(stan.data, df$const[!duplicated(df$adm2)], betas, gammas, rep(1, stan.data$M)))
    print(calc.likeli(stan.data, df$const[!duplicated(df$adm2)], betas + rnorm(4), gammas, rep(1, stan.data$M)))
    print(calc.likeli(stan.data, df$const[!duplicated(df$adm2)], betas, gammas + rnorm(1), rep(1, stan.data$M)))

    print(calc.likeli.partial(stan.data, betas, gammas, rep(1, stan.data$M)))
    print(calc.likeli.partial(stan.data, betas + rnorm(4), gammas, rep(1, stan.data$M)))
    print(calc.likeli.partial(stan.data, betas, gammas + rnorm(1), rep(1, stan.data$M)))
}
