##setwd("~/research/gcp/hierarchical-estimation/logspec-ols")

do.checks <- F

calc.expected.demeaned <- function(K, dmxx, dmyy, zz, mm, betas, gammas) {
    obsmean <- 0
    for (kk in 1:K)
        obsmean <- obsmean + betas[kk] * dmxx[kk, ] * exp(zz[kk, , ] %*% gammas[kk, ])

    obsmean
}

calc.likeli.demeaned <- function(K, dmxx, dmyy, zz, mm, betas, gammas, sigma) {
    obsmean <- calc.expected.demeaned(K, dmxx, dmyy, zz, mm, betas, gammas)
    sum(dnorm(dmyy, obsmean, sigma[mm], log=T))
}

source("../example/logspec-data.R", chdir=T)

df$dmrate <- regional.demean(df$rate, df$adm2)
## Demean all predictors
df$dmbin1 <- regional.demean(df$bin1, df$adm2)
df$dmbin2 <- regional.demean(df$bin2, df$adm2)
df$dmbin4 <- regional.demean(df$bin4, df$adm2)
df$dmbin5 <- regional.demean(df$bin5, df$adm2)

if (do.checks) {
    ## Check that true values are optimum
    teff0 <- .1 # Size of temperature effect
    tstar <- 21 # Lowest-effect temperature
    betas <- teff0 * (c(12, 16.5, 25.5, 30) - tstar)^2 - teff0 * (21 - tstar)^2
    gammas <- t(matrix(rep(c(-.005, -.05), 4), 2, 4))

    df$dmexprate <- betas[1] * exp(-df$log_gdppc / 20 - df$meant / 200) * df$dmbin1 + betas[2] * exp(-df$log_gdppc / 20 - df$meant / 200) * df$dmbin2 + betas[3] * exp(-df$log_gdppc / 20 - df$meant / 200) * df$dmbin4 + betas[4] * exp(-df$log_gdppc / 20 - df$meant / 200) * df$dmbin5 # Exact match
}

K <- stan.data$K
L <- stan.data$L
dmxx <- t(as.matrix(df[, c('dmbin1', 'dmbin2', 'dmbin4', 'dmbin5')]))
dmyy <- df$dmrate
zz <- stan.data$zz

if (do.checks) {
    df$dmexprate <- calc.expected.demeaned(K, dmxx, dmyy, zz, mm, betas, gammas)
}

## Approximation 1: No covariate effect
stacked <- lm(dmrate ~ 0 + dmbin1 + dmbin2 + dmbin4 + dmbin5, data=df)

## Approximation 2: Identically distributed
objective <- function(params) {
    betas <- params[1:K]
    gammas <- matrix(params[(K+1):((L+1)*K)], K, L)
    sigma <- rep(params[((L+1)*K+1)], stan.data$M)
    -calc.likeli.demeaned(K, dmxx, dmyy, zz, mm, betas, gammas, sigma)
}

params <- c(stacked$coeff, rep(0, K*L), sd(stacked$residuals))
soln <- optim(params, objective)

## Approximation 3: State-clustered errors
gammas <- matrix(params[(K+1):((L+1)*K)], K, L)
df$dmbin1exp <- df$dmbin1 * exp(gammas[1,1] * df$meant + gammas[1,2] * df$log_gdppc)
df$dmbin2exp <- df$dmbin2 * exp(gammas[2,1] * df$meant + gammas[2,2] * df$log_gdppc)
df$dmbin4exp <- df$dmbin4 * exp(gammas[3,1] * df$meant + gammas[3,2] * df$log_gdppc)
df$dmbin5exp <- df$dmbin5 * exp(gammas[4,1] * df$meant + gammas[4,2] * df$log_gdppc)
stacked <- lm(dmrate ~ 0 + dmbin1exp + dmbin2exp + dmbin4exp + dmbin5exp, data=df)

sigma <- c()
for (jj in 1:stan.data$M)
    sigma <- c(sigma, sd(stacked$residuals[df$adm1 == jj]))

objective2 <- function(params) {
    betas <- params[1:K]
    gammas <- matrix(params[(K+1):((L+1)*K)], K, L)
    sigma <- params[((L+1)*K+1):length(params)]
    -calc.likeli.demeaned(K, dmxx, dmyy, zz, mm, betas, gammas, sigma)
}

params <- c(stacked$coeff, soln$par[(K+1):((L+1)*K)], sigma)
soln <- optim(params, objective2)

## Get Hessian by by assuming sigma
sigma <- soln$par[((L+1)*K+1):length(soln$par)]

objective3 <- function(params) {
    betas <- params[1:K]
    gammas <- matrix(params[(K+1):((L+1)*K)], K, L)
    -calc.likeli.demeaned(K, dmxx, dmyy, zz, mm, betas, gammas, sigma)
}

params <- c(stacked$coeff, rep(0, K*L))
soln <- optim(params, objective3, hessian=T)

betas <- soln$par[1:K]
gammas <- soln$par[(K+1):((L+1)*K)]

ses <- sqrt(abs(diag(solve(soln$hessian))))
