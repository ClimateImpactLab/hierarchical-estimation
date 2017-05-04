##setwd("~/research/gcp/hierarchical-estimation/ranges")

source("../interpolate/logspec.R")

df <- read.csv("../example/true-binned.csv")
df$log_gdppc <- log(df$gdppc)

## Find the coefficients
result <- estimate.logspec(df$rate, df[, c('bin1', 'bin2', 'bin4', 'bin5')],
                           df[!duplicated(df$adm1), c('meant', 'log_gdppc', 'meant', 'log_gdppc',
                                                      'meant', 'log_gdppc', 'meant', 'log_gdppc')],
                           df$adm1, df$adm2)

calc.likeli.nosigma <- function(K, L, dmxxs, dmyy, zzs, mm, betas, gammas) {
    dmyy.exp <- calc.expected.demeaned(K, L, dmxxs, dmyy, zzs, mm, betas, gammas)

    sigmas <- c()
    for (jj in unique(mm)) {
        included <- mm == jj
        residuals <- dmyy.exp[included] - dmyy[included]
        sigmas <- c(sigmas, sd(residuals))
    }

    calc.likeli.demeaned(K, L, dmxxs, dmyy, zzs, mm, betas, gammas, sigmas)
}

dmyy <- regional.demean(df$rate, df$adm2)
dmxxs <- df[, c('bin1', 'bin2', 'bin4', 'bin5')]
for (kk in 1:4)
    dmxxs[, kk] <- regional.demean(dmxxs[, kk], df$adm2)

zzs <- df[!duplicated(df$adm1), c('meant', 'log_gdppc', 'meant', 'log_gdppc',
                                  'meant', 'log_gdppc', 'meant', 'log_gdppc')]

calc.likeli.nosigma(4, 2, dmxxs, dmyy, zzs, df$adm1, result$betas, result$gammas)

## Use Metropolis-Hastings with 1 seeds
## accept.target <- .44 # Optimal for 1-D sampling
## Want (a^44) * (b^56) = 1; Say a = 1.1, (1 / (1.1^44))^(1 / 56) = 0.9278487

methast <- function(K, L, dmxxs, dmyy, zzs, mm, iter, beta0, gamma0, betaerr, gammaerr) {
    params = matrix(NA, iter, K + K*L)
    params[1, ] = c(beta0, gamma0)

    last.likeli <- calc.likeli.nosigma(K, L, dmxxs, dmyy, zzs, mm, beta0, gamma0)
    sd.product <- 1

    for (ii in 2:iter) {
        if (ii %% 100 == 0)
            print(ii)
        beta.sample <- rnorm(K, params[ii-1, 1:K], betaerr * sd.product)
        gamma.sample <- matrix(rnorm(K*L, params[ii-1, (K+1):(K+K*L)], c(gammaerr) * sd.product), K, L)

        this.likeli <- calc.likeli.nosigma(K, L, dmxxs, dmyy, zzs, mm, beta.sample, gamma.sample)
        prob <- exp(this.likeli - last.likeli)

        if (min(prob, 1) > runif(1)) {
            params[ii, ] <- c(beta.sample, gamma.sample)
            last.likeli <- this.likeli
            sd.product <- sd.product * 1.1 # was too modest
        } else {
            params[ii, ] <- params[ii-1, ]
            sd.product <- sd.product * 0.9278487 # was too bold
        }
    }

    return(params)
}

params <- methast(4, 2, dmxxs, dmyy, zzs, df$adm1, 1000, result$betas, result$gammas, result$ses.betas, result$ses.gammas)

cov(params)
apply(params, 2, function(xx) quantile(xx))

params <- methast(4, 2, dmxxs, dmyy, zzs, df$adm1, 1000, result$betas, result$gammas, rep(1, 4), matrix(1, 4, 2))
apply(params, 2, sd)

## Use Metropolis-Hastings with N seeds
repeated.methast <- function(K, L, dmxxs, dmyy, zzs, mm, iter, seeds, beta0, gamma0, betaerr, gammaerr) {
    params = matrix(NA, 0, K + K*L)

    for (seed in 1:seeds) {
        print(paste("Seed", seed))
        params.seed <- methast(4, 2, dmxxs, dmyy, zzs, df$adm1, iter, result$betas, result$gammas, result$ses.betas, result$ses.gammas)
        params <- rbind(params, params.seed)
    }

    params
}

params <- repeated.methast(4, 2, dmxxs, dmyy, zzs, df$adm1, 1000, 4, result$betas, result$gammas, result$ses.betas, result$ses.gammas)

cov(params)
apply(params, 2, function(xx) quantile(xx, probs=c(.025, .1, .5, .9, .975)))
colMeans(params)

## Also, calculate effective SD that would fit the 95% range
sd.ols <- c(result$ses.betas, result$ses.gammas)
sd.bayes <- apply(params, 2, sd)

obsd <- apply(apply(params, 2, function(xx) quantile(xx, probs=c(.025, .975))), 2, diff)
expd <- qnorm(.975) - qnorm(.025)
sd.tail <- obsd / expd

sd.conservative <- pmax(sd.ols, sd.bayes, sd.tail)
