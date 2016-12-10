library(nnls)
library(RcppArmadillo)

calc.expected.demeaned <- function(K, L, dmxxs, dmyy, zzs, mm, betas, gammas) {
    obsmean <- 0
    for (kk in 1:K)
        obsmean <- obsmean + betas[kk] * dmxxs[, kk] * exp(as.matrix(zzs[, ((kk-1)*L + 1):(kk*L)]) %*% gammas[kk, ])[mm]

    obsmean
}

calc.likeli.demeaned <- function(K, L, dmxxs, dmyy, zzs, mm, betas, gammas, sigma) {
    obsmean <- calc.expected.demeaned(K, L, dmxxs, dmyy, zzs, mm, betas, gammas)
    sum(dnorm(dmyy, obsmean, sigma[mm], log=T))
}

regional.demean <- function(values, regions) {
    for (region in unique(regions)) {
        regioniis <- which(regions == region)
        values[regioniis] <- values[regioniis] - mean(values[regioniis])
    }

    values
}

check.arguments <- function(yy, xxs, zzs, adm1, adm2) {
    if (is.null(nrow(xxs)) || is.null(ncol(xxs)))
        stop("xxs must be a matrix or data.frame.")

    if (is.null(nrow(zzs)) || is.null(ncol(zzs)))
        stop("zzs must be a matrix or data.frame.")

    N <- length(yy)
    if (nrow(xxs) != N || length(adm1) != N || length(adm2) != N)
        stop("yy, xxs, adm1, and adm2 must all have the same number of observations.")

    K <- ncol(xxs)
    L <- ncol(zzs) / K

    if (L %% 1 != 0)
        stop("zzs must have columns a whole number multiple of xxs's")

    M <- max(adm1)
    if (length(unique(adm1)) != M)
        stop("adm1 must contain intgers from 1 to M.")

    if (nrow(zzs) != M)
        stop("zzs must have the same number of rows as ADM1 values.")

    return(list(K=K, L=L, N=N, M=M))
}

demean.yxs <- function(yy, xxs, adm2) {
    ## De-mean observations
    dmyy <- regional.demean(yy, adm2)
    dmxxs <- xxs
    for (kk in 1:K)
        dmxxs[, kk] <- regional.demean(dmxxs[, kk], adm2)

    list(dmyy=dmyy, dmxxs=dmxxs)
}

estimate.logspec <- function(yy, xxs, zzs, adm1, adm2, maxiter=1000, initgammas=NULL) {
    list2env(check.arguments(yy, xxs, zzs, adm1, adm2), parent.frame())
    list2env(demean.yxs(yy, xxs, adm2), parent.frame())
    dmxxs.orig <- dmxxs

    print("Iterating...")

    bestlikeli <- -Inf # best likelihood we've seen
    bestgammas <- matrix(0, K, L) # gammas corresponding to bestlikeli
    armijo.factor <- 1 # the amount of movement away from bestgammas

    ## Start with no covariate effects
    dmxxs <- dmxxs.orig

    if (is.null(initgammas)) {
        gammas <- bestgammas
    } else {
        gammas <- initgammas
        for (kk in 1:K)
            dmxxs[, kk] <- dmxxs[, kk] * exp(as.matrix(zzs[, ((kk-1)*L + 1):(kk*L)]) %*% gammas[kk, ])[adm1]
    }

    ## Weighted by region-error by dividing by standard error
    dmyy.weighted <- dmyy
    dmxxs.weighted <- dmxxs

    for (iter in 1:maxiter) {
        ## Perform stacked regression to get signs
        stacked <- fastLm(dmxxs.weighted, dmyy.weighted)
        betas <- stacked$coeff

        ## Check the performance of this stacked set of betas
        stage1.sigma <- rep(NA, M) # Use same variable as stage1-2, in case we exit here
        for (jj in 1:M) {
            included <- adm1 == jj
            stage1.sigma[jj] <- mean(sd(stacked$residuals) * dmyy[included] / dmyy.weighted[included])
        }

        likeli1 <- calc.likeli.demeaned(K, L, dmxxs.orig, dmyy, zzs, adm1, betas, gammas, stage1.sigma)

        ## Check if we have converged
        if (abs(likeli1 - bestlikeli) < 1e-6 || armijo.factor < 1e-6)
            break

        if (likeli1 > bestlikeli) {
            bestlikeli <- likeli1
            bestgammas <- gammas
        }

        ## Perform a regression on each state
        stage1.betas <- matrix(NA, M, K)
        stage1.sigma <- rep(NA, M)
        for (jj in 1:M) {
            included <- adm1 == jj
            dmxxsjj <- as.matrix(dmxxs[included,])
            modjj <- nnnpls(dmxxsjj, dmyy[included], stacked$coeff)

            ## Multiply back in exponent, if dmbins were generated with an assumed gamma
            for (kk in 1:K)
                stage1.betas[jj, kk] <- modjj$x[kk] * mean(dmxxsjj[, kk] / dmxxs.orig[included, kk], na.rm=T)

            stage1.sigma[jj] <- sd(modjj$residuals) # NOTE: residual standard error is not quite this
        }

        ## Prepare values for second stage regressions
        stage1.logbetas <- log(abs(stage1.betas))
        stage1.logbetas[stage1.logbetas == -Inf] <- NA

        ## Second stage regression for each predictor
        betas <- rep(NA, K)
        gammas <- matrix(NA, K, L)
        for (kk in 1:K) {
            valid <- is.finite(stage1.logbetas[, kk])

            modkk <- fastLm(cbind(rep(1, sum(valid)), zzs[valid, ((kk-1)*L + 1):(kk*L)]), stage1.logbetas[valid, kk])

            ## Prepare values for log likelihood
            betas[kk] <- exp(modkk$coeff[1]) * sign(stacked$coeff[kk])
            gammas[kk, ] <- modkk$coeff[-1]
        }

        likeli2 <- calc.likeli.demeaned(K, L, dmxxs.orig, dmyy, zzs, adm1, betas, gammas, stage1.sigma)

        ## Report progress
        print(c(iter, log2(1/armijo.factor), as.integer(max(likeli1, likeli2))))

        ## Check if we have converged
        if (abs(likeli2 - bestlikeli) < 1e-6 || armijo.factor < 1e-6)
            break

        if (likeli2 > bestlikeli) {
            bestlikeli <- likeli2
            bestgammas <- gammas
        } else {
            ## Oops!  We spun too far out
            gammas <- (bestgammas + armijo.factor * gammas) / (1 + armijo.factor)
            armijo.factor <- armijo.factor / 2
            ## Note: We may never achive as good a likelihood, since it may be inconsistent with the sigmas from stage1
        }

        dmxxs <- dmxxs.orig
        for (kk in 1:K)
            dmxxs[, kk] <- dmxxs[, kk] * exp(as.matrix(zzs[, ((kk-1)*L + 1):(kk*L)]) %*% gammas[kk, ])[adm1]
        for (jj in 1:M) {
            included <- adm1 == jj
            dmyy.weighted[included] <- dmyy[included] / stage1.sigma[jj]
            dmxxs.weighted[included,] <- dmxxs[included,] / stage1.sigma[jj]
        }
    }

    list(betas=betas, gammas=gammas, sigmas=stage1.sigma)
}

stacked.regression <- function(L, gammas, dmyy, dmxxs.orig, zzs, mm) {
    dmxxs <- dmxxs.orig
    for (kk in 1:K)
        dmxxs[, kk] <- dmxxs[, kk] * exp(as.matrix(zzs[, ((kk-1)*L + 1):(kk*L)]) %*% gammas[kk, ])[mm]

    fastLm(dmxxs, dmyy)
}

## Get stacked betas
stacked.betas <- function(L, gammas, dmyy, dmxxs.orig, zzs, mm) {
    stacked <- stacked.regression(L, gammas, dmyy, dmxxs.orig, zzs, mm)

    stacked$coeff
}

estimate.logspec.optim <- function(yy, xxs, zzs, adm1, adm2, initgammas=NULL) {
    list2env(check.arguments(yy, xxs, zzs, adm1, adm2), parent.frame())
    list2env(demean.yxs(yy, xxs, adm2), parent.frame())
    dmxxs.orig <- dmxxs

    if (is.null(initgammas)) {
        ## Approximation 1: No covariate effect
        stacked <- fastLm(dmxxs.orig, dmyy)
        sigma <- sd(stacked$residuals)

        betas <- stacked$coeff
        gammas <- matrix(0, K, L)
    } else {
        gammas <- initgammas
        stacked <- stacked.regression(L, gammas, dmyy, dmxxs.orig, zzs, adm1)
        betas <- stacked$coeff
        sigma <- sd(stacked$residuals)
    }

    print(c("Step 1:", calc.likeli.demeaned(K, L, dmxxs.orig, dmyy, zzs, adm1, betas, gammas, rep(sigma, M))))

    ## Approximation 2: Identically distributed
    objective <- function(params) {
        gammas <- matrix(params[1:(L*K)], K, L)
        sigma <- rep(params[L*K+1], M)

        betas <- stacked.betas(L, gammas, dmyy, dmxxs.orig, zzs, adm1)

        -calc.likeli.demeaned(K, L, dmxxs.orig, dmyy, zzs, adm1, betas, gammas, sigma)
    }

    params <- c(gammas, sigma)
    soln <- optim(params, objective)

    gammas <- matrix(soln$par[1:(L*K)], K, L)
    sigma <- soln$par[L*K+1]

    betas <- stacked.betas(L, gammas, dmyy, dmxxs.orig, zzs, adm1)

    print(c("Step 2:", calc.likeli.demeaned(K, L, dmxxs.orig, dmyy, zzs, adm1, betas, gammas, rep(sigma, M))))

    ## Approximation 3: State-clustered errors
    dmxxs <- dmxxs.orig
    for (kk in 1:K)
        dmxxs[, kk] <- dmxxs[, kk] * exp(as.matrix(zzs[, ((kk-1)*L + 1):(kk*L)]) %*% gammas[kk, ])[adm1]
    stacked <- fastLm(dmxxs, dmyy)

    betas <- stacked$coeff

    sigma <- c()
    for (jj in 1:M)
        sigma <- c(sigma, sd(stacked$residuals[adm1 == jj]))

    print(c("Step 3:", calc.likeli.demeaned(K, L, dmxxs.orig, dmyy, zzs, adm1, betas, gammas, sigma)))

    objective2 <- function(params) {
        gammas <- matrix(params[1:(L*K)], K, L)
        sigma <- params[(L*K+1):length(params)]

        betas <- stacked.betas(L, gammas, dmyy, dmxxs.orig, zzs, adm1)

        -calc.likeli.demeaned(K, L, dmxxs.orig, dmyy, zzs, adm1, betas, gammas, sigma)
    }

    params <- c(soln$par[1:(L*K)], sigma)
    soln <- optim(params, objective2)

    gammas <- matrix(soln$par[1:(L*K)], K, L)
    sigma <- soln$par[(L*K+1):length(soln$par)]

    betas <- stacked.betas(L, gammas, dmyy, dmxxs.orig, zzs, adm1)

    print(c("Step 4:", calc.likeli.demeaned(K, L, dmxxs.orig, dmyy, zzs, adm1, betas, gammas, sigma)))

    list(betas=soln$par[1:K], gammas=matrix(soln$par[(K+1):((L+1)*K)], K, L), ses.betas=ses[1:K], ses.gammas=matrix(ses[(K+1):((L+1)*K)], K, L), sigma=sigma)
}

source("methast.R")
source("ranges.R")
