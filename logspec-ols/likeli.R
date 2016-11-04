##setwd("~/projects/gcp/hierarchical-estimation/logspec-ols")

do.tests <- F

calc.likeli <- function(stan.data, alphas, betas, gammas, sigma) {
    obsmean <- alphas[stan.data$nn]
    for (kk in 1:stan.data$K)
        obsmean <- obsmean + betas[kk] * stan.data$xx[kk, ] * exp(stan.data$zz[kk, , ] %*% gammas[kk, ])

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
    tstar <- 22 # Lowest-effect temperature
    betas <- teff0 * (c(10.15975, 17.64027, 24.86336, 31.78611) - tstar)^2
    gammas <- t(matrix(rep(c(-.05, -.005), 4), 2, 4))

    print(calc.likeli(stan.data, df$const[!duplicated(df$adm2)], betas, gammas, rep(1, stan.data$M)))
    print(calc.likeli(stan.data, df$const[!duplicated(df$adm2)], betas + rnorm(4), gammas, rep(1, stan.data$M)))
    print(calc.likeli(stan.data, df$const[!duplicated(df$adm2)], betas, gammas + rnorm(1), rep(1, stan.data$M)))

    print(calc.likeli.partial(stan.data, betas, gammas, rep(1, stan.data$M)))
    print(calc.likeli.partial(stan.data, betas + rnorm(4), gammas, rep(1, stan.data$M)))
    print(calc.likeli.partial(stan.data, betas, gammas + rnorm(1), rep(1, stan.data$M)))
}
