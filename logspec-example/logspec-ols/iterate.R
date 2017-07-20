##setwd("~/research/gcp/hierarchical-estimation/logspec-ols")

library(nnls)

do.checks <- F

calc.expected.demeaned <- function(K, dmxx, zz, mm, betas, gammas) {
    obsmean <- 0
    for (kk in 1:K)
        obsmean <- obsmean + betas[kk] * dmxx[kk, ] * exp(zz[kk, , ] %*% gammas[kk, ])

    obsmean
}

calc.likeli.demeaned <- function(K, dmxx, dmyy, zz, mm, betas, gammas, sigma) {
    obsmean <- calc.expected.demeaned(K, dmxx, zz, mm, betas, gammas)
    sum(dnorm(dmyy, obsmean, sigma[mm], log=T))
}

source("../example/logspec-data.R", chdir=T)

df$dmrate <- regional.demean(df$rate, df$adm2)
## Demean all predictors
df$dmbin1.orig <- regional.demean(df$bin1, df$adm2)
df$dmbin2.orig <- regional.demean(df$bin2, df$adm2)
df$dmbin4.orig <- regional.demean(df$bin4, df$adm2)
df$dmbin5.orig <- regional.demean(df$bin5, df$adm2)

## Make copies of bin values, for adjusting
df$dmbin1 <- df$dmbin1.orig
df$dmbin2 <- df$dmbin2.orig
df$dmbin4 <- df$dmbin4.orig
df$dmbin5 <- df$dmbin5.orig

bestlikeli <- -Inf
bestgammas <- rep(0, stan.data$K*stan.data$L)
armijo.factor <- 1
for (iter in 1:1000) {
    ## Perform stacked regression to get signs
    stacked <- lm(dmrate ~ 0 + dmbin1 + dmbin2 + dmbin4 + dmbin5, data=df)

    ## Perform a regression on each state
    stage1 <- data.frame(adm1=c(), beta1=c(), beta2=c(), beta3=c(), beta4=c(), meant=c(), log_gdppc=c(), sigma=c())
    for (jj in unique(df$adm1)) {
        subdf <- subset(df, adm1 == jj)
        modjj <- nnnpls(as.matrix(subdf[, c('dmbin1', 'dmbin2', 'dmbin4', 'dmbin5')]), subdf$dmrate, stacked$coeff)
        ## Multiply back in exponent, if dmbins were generated with an assumed gamma
        beta1 <- modjj$x[1] * mean(subdf$dmbin1 / subdf$dmbin1.orig, na.rm=T)
        beta2 <- modjj$x[2] * mean(subdf$dmbin2 / subdf$dmbin2.orig, na.rm=T)
        beta4 <- modjj$x[3] * mean(subdf$dmbin4 / subdf$dmbin4.orig, na.rm=T)
        beta5 <- modjj$x[4] * mean(subdf$dmbin5 / subdf$dmbin5.orig, na.rm=T)

        stage1 <- rbind(stage1, data.frame(adm1=jj, beta1, beta2, beta4, beta5, meant=subdf$meant[1], log_gdppc=subdf$log_gdppc[1], sigma=sd(modjj$residuals))) # NOTE: residual standard error is not quite this
    }

    ## Prepare values for second stage regressions
    stage1$logbeta1 <- log(abs(stage1$beta1))
    stage1$logbeta2 <- log(abs(stage1$beta2))
    stage1$logbeta4 <- log(abs(stage1$beta4))
    stage1$logbeta5 <- log(abs(stage1$beta5))

    stage1$logbeta1[stage1$logbeta1 == -Inf] = NA
    stage1$logbeta2[stage1$logbeta2 == -Inf] = NA
    stage1$logbeta4[stage1$logbeta4 == -Inf] = NA
    stage1$logbeta5[stage1$logbeta5 == -Inf] = NA

    ## Second stage regression for each predictor
    mod1 <- lm(logbeta1 ~ meant + log_gdppc, data=stage1)
    mod2 <- lm(logbeta2 ~ meant + log_gdppc, data=stage1)
    mod4 <- lm(logbeta4 ~ meant + log_gdppc, data=stage1)
    mod5 <- lm(logbeta5 ~ meant + log_gdppc, data=stage1)

    ## Prepare values for log likelihood
    betas <- exp(c(mod1$coeff[1], mod2$coeff[1], mod4$coeff[1], mod5$coeff[1])) * sign(stacked$coeff)
    gammas <- t(matrix(c(mod1$coeff[-1], mod2$coeff[-1], mod4$coeff[-1], mod5$coeff[-1]), 2, 4))

    dmxx <- t(as.matrix(df[, c('dmbin1.orig', 'dmbin2.orig', 'dmbin4.orig', 'dmbin5.orig')]))
    likeli <- calc.likeli.demeaned(stan.data$K, dmxx, df$dmrate, stan.data$zz, stan.data$mm, betas, gammas, stage1$sigma)
    print(c(iter, likeli))

    if (abs(likeli - bestlikeli) < 1e-6 || armijo.factor < 1e-6)
        break

    if (likeli > bestlikeli) {
        bestlikeli <- likeli
        bestgammas <- gammas
    } else {
        ## Oops!  We spun too far out
        gammas <- (bestgammas + armijo.factor * gammas) / (1 + armijo.factor)
        armijo.factor <- armijo.factor / 2
        ## Note: We may never achive as good a likelihood, since it may be inconsistent with the sigmas from stage1
    }

    df$dmbin1 <- df$dmbin1.orig * exp(stan.data$zz[1, , ] %*% gammas[1, ])
    df$dmbin2 <- df$dmbin2.orig * exp(stan.data$zz[2, , ] %*% gammas[2, ])
    df$dmbin4 <- df$dmbin4.orig * exp(stan.data$zz[3, , ] %*% gammas[3, ])
    df$dmbin5 <- df$dmbin5.orig * exp(stan.data$zz[4, , ] %*% gammas[4, ])
}
