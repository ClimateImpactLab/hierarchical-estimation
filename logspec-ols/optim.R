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
## Demean all predictors
df$dmbin1 <- regional.demean(df$bin1, df$adm2)
df$dmbin2 <- regional.demean(df$bin2, df$adm2)
df$dmbin4 <- regional.demean(df$bin4, df$adm2)
df$dmbin5 <- regional.demean(df$bin5, df$adm2)

K <- stan.data$K
L <- stan.data$L
dmxx <- as.matrix(df[, c('dmbin1', 'dmbin2', 'dmbin4', 'dmbin5')])
dmyy <- df$dmrate
zz <- stan.data$zz

objective <- function(params) {
    betas <- params[1:K]
    gammas <- matrix(params[(K+1):((L+1)*K)], K, L)
    sigma <- params[((L+1)*K+1):length(params)]
    -calc.likeli.demeaned(K, dmxx, dmyy, zz, mm, betas, gammas, sigma)
}

stacked <- lm(dmrate ~ 0 + dmbin1 + dmbin2 + dmbin4 + dmbin5, data=df)

params <- c(stacked$coeff, rep(0, K*L), rep(1, stan.data$M))
soln <- optim(params, objective)

sigma <- soln$par[((L+1)*K+1):length(soln$par)]

objective2 <- function(params) {
    betas <- params[1:K]
    gammas <- matrix(params[(K+1):((L+1)*K)], K, L)
    -calc.likeli.demeaned(K, dmxx, dmyy, zz, mm, betas, gammas, sigma)
}

soln2 <- optim(soln$par[1:((L+1)*K)], objective2, hessian=T)

betas2 <- soln2$par[1:K]
gammas2 <- soln2$par[(K+1):((L+1)*K)]
ses <- sqrt(abs(diag(solve(soln2$hessian))))
