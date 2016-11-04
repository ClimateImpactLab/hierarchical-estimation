##setwd("~/research/gcp/hierarchical-estimation/logspec-ols")

library(nnls)

calc.likeli.demeaned <- function(K, dmxx, dmyy, zz, mm, betas, gammas, sigma) {
    obsmean <- 0
    for (kk in 1:K)
        obsmean <- obsmean + betas[kk] * dmxx[kk, ] * exp(zz[kk, , ] %*% gammas[kk, ])

    sum(dnorm(dmyy, obsmean, sigma[mm], log=T))
}

source("../example/logspec-data.R", chdir=T)

df$dmrate <- regional.demean(df$rate, df$adm2)
## Make copies of bin values, for adjusting
df$xxbin1 <- df$bin1
df$xxbin2 <- df$bin2
df$xxbin4 <- df$bin4
df$xxbin5 <- df$bin5

prevlikeli <- -Inf
for (iter in 1:1000) {
    ## Demean all predictors
    df$dmbin1 <- regional.demean(df$xxbin1, df$adm2)
    df$dmbin2 <- regional.demean(df$xxbin2, df$adm2)
    df$dmbin4 <- regional.demean(df$xxbin4, df$adm2)
    df$dmbin5 <- regional.demean(df$xxbin5, df$adm2)

    ## Perform stacked regression to get signs
    stacked <- lm(dmrate ~ 0 + dmbin1 + dmbin2 + dmbin4 + dmbin5, data=df)

    ## Perform a regression on each state
    stage1 <- data.frame(adm1=c(), beta1=c(), beta2=c(), beta3=c(), beta4=c(), meant=c(), log_gdppc=c(), sigma=c())
    for (jj in unique(df$adm1)) {
        subdf <- subset(df, adm1 == jj)
        modjj <- nnnpls(as.matrix(subdf[, c('dmbin1', 'dmbin2', 'dmbin4', 'dmbin5')]), subdf$dmrate, stacked$coeff)
        stage1 <- rbind(stage1, data.frame(adm1=jj, beta1=modjj$x[1], beta2=modjj$x[2], beta4=modjj$x[3], beta5=modjj$x[4], meant=subdf$meant[1], log_gdppc=subdf$log_gdppc[1], sigma=sd(modjj$residuals))) # NOTE: residual standard error is not quite this
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

    dmxx <- as.matrix(df[, c('dmbin1', 'dmbin2', 'dmbin4', 'dmbin5')])
    likeli <- calc.likeli.demeaned(stan.data$K, dmxx, df$dmrate, stan.data$zz, stan.data$mm, betas, gammas, stage1$sigma)
    print(c(iter, likeli))

    if (abs(likeli - prevlikeli) < 1e-6)
        break

    prevlikeli <- likeli
    df$xxbin1 <- df$bin1 * exp(stan.data$zz[1, , ] %*% gammas[1, ])
    df$xxbin2 <- df$bin2 * exp(stan.data$zz[2, , ] %*% gammas[2, ])
    df$xxbin4 <- df$bin4 * exp(stan.data$zz[3, , ] %*% gammas[3, ])
    df$xxbin5 <- df$bin5 * exp(stan.data$zz[4, , ] %*% gammas[4, ])
}
