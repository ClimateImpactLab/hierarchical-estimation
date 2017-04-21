library(nnls)
library(RcppArmadillo)

## Calculate x_k exp(sum gamma_kl z_l) as a NxK matrix
calc.covariated.predictors <- function(dmxxs, zzs, kls, mm, gammas) {
    result <- dmxxs # Only modify this for predictors with covariates
    gammas.so.far <- 0 # Keep track of how many coefficients used

    for (kk in 1:nrow(kls)) {
        gammas.here <- sum(kls[kk, ])
        if (gammas.here == 0)
            next # Nothing to do: already dmxxs[, kk]

        mygammas <- gammas[(gammas.so.far+1):(gammas.so.far+gammas.here)]
        gammas.so.far <- gammas.so.far + gammas.here
        result[, kk] <- dmxxs[, kk] * exp(as.matrix(zzs[, kls[kk, ]]) %*% mygammas)[mm]
    }

    result
}

## Calculate the expected value, after it has been demeaned
calc.expected.demeaned <- function(dmxxs, dmyy, zzs, kls, mm, betas, gammas) {
    covariated.predictors <- calc.covariated.predictors(dmxxs, zzs, kls, mm, gammas)

    obsmean <- 0
    for (kk in 1:nrow(kls))
        obsmean <- obsmean + betas[kk] * covariated.predictors[, kk]

    obsmean
}

## Calculate the likelihood of the given parameters, after observations have been demeaned
calc.likeli.demeaned <- function(dmxxs, dmyy, zzs, kls, mm, betas, gammas, sigma, weights) {
    obsmean <- calc.expected.demeaned(dmxxs, dmyy, zzs, kls, mm, betas, gammas)
    sum(weights * dnorm(dmyy, obsmean, sigma[mm], log=T))
}

## Demean a set of values by region (partial out region fixed effects)
regional.demean <- function(values, regions) {
    for (region in unique(regions)) {
        regioniis <- which(regions == region)
        values[regioniis] <- values[regioniis] - mean(values[regioniis])
    }

    values
}

## Check that all of the given arguments are correctly specified
check.arguments <- function(yy, xxs, zzs, kls, adm1, adm2) {
    if (is.null(nrow(xxs)) || is.null(ncol(xxs)))
        stop("xxs must be a matrix or data.frame.")

    if (is.null(nrow(zzs)) || is.null(ncol(zzs)))
        stop("zzs must be a matrix or data.frame.")

    N <- length(yy)
    if (nrow(xxs) != N || length(adm1) != N || length(adm2) != N)
        stop("yy, xxs, adm1, and adm2 must all have the same number of observations.")

    K <- ncol(xxs)
    L <- ncol(zzs)

    if (nrow(kls) != K)
        stop("kls must have as many rows as xxs has columns.")

    if (ncol(kls) != L)
        stop("kls must have as many columns as zzs.")

    if (!is.logical(kls))
        stop("kls must consist only of trues and falses.")

    M <- max(adm1)
    if (length(unique(adm1)) != M)
        stop("adm1 must contain intgers from 1 to M.")

    if (nrow(zzs) != M)
        stop("zzs must have the same number of rows as ADM1 values.")

    return(list(K=K, L=L, N=N, M=M))
}

## Demean all observations by ADM2 regions (partial out ADM2 fixed effects)
demean.yxs <- function(K, yy, xxs, adm2) {
    ## De-mean observations
    dmyy <- regional.demean(yy, adm2)
    dmxxs <- xxs
    for (kk in 1:K)
        dmxxs[, kk] <- regional.demean(dmxxs[, kk], adm2)

    list(dmyy=dmyy, dmxxs=dmxxs)
}

## Estimate the parameters of the log specification
estimate.logspec <- function(yy, xxs, zzs, kls, adm1, adm2, weights=1, maxiter=1000, initgammas=NULL) {
    list2env(check.arguments(yy, xxs, zzs, kls, adm1, adm2), environment())
    list2env(demean.yxs(K, yy, xxs, adm2), environment())
    dmxxs.orig <- dmxxs

    print("Iterating...")

    bestlikeli <- -Inf # best likelihood we've seen
    bestgammas <- rep(0, sum(kls)) # gammas corresponding to bestlikeli
    bestbetas <- rep(0, K) # betas corresponding to bestlikeli
    bestsigmas <- rep(Inf, M) # adm1.sigma corresponding to bestlikeli
    armijo.factor <- 1 # the amount of movement away from bestgammas
    bestgravity <- F # should we move toward best this iteration

    ## Start with no covariate effects
    dmxxs <- dmxxs.orig

    if (is.null(initgammas)) {
        gammas <- bestgammas
    } else {
        gammas <- initgammas
        dmxxs <- calc.covariated.predictors(dmxxs, zzs, kls, adm1, gammas)
    }

    ## Weighted by region-error by dividing by standard error
    dmyy.weighted <- dmyy * sqrt(weights)
    dmxxs.weighted <- dmxxs
    for (kk in 1:K)
        dmxxs.weighted[, kk] <- dmxxs[, kk] * sqrt(weights)

    for (iter in 1:maxiter) {
        ## Perform stacked regression to get signs
        stacked <- fastLm(dmxxs.weighted, dmyy.weighted)
        betas <- stacked$coeff

        ## Check the performance of this stacked set of betas
        adm1.sigma <- rep(NA, M)
        for (jj in 1:M) {
            included <- adm1 == jj
            adm1.sigma[jj] <- mean(sd(stacked$residuals[included]) * dmyy[included] / dmyy.weighted[included])
        }

        if (bestgravity) {
            betas <- (bestbetas + armijo.factor * betas) / (1 + armijo.factor)
            adm1.sigma <- (bestsigmas + armijo.factor * adm1.sigma) / (1 + armijo.factor)
        }

        likeli1 <- calc.likeli.demeaned(dmxxs.orig, dmyy, zzs, kls, adm1, betas, gammas, adm1.sigma, weights)

        ## Check if we have converged
        if (abs(likeli1 - bestlikeli) < 1e-6 || armijo.factor < 1e-6)
            break

        if (likeli1 > bestlikeli) {
            bestlikeli <- likeli1
            bestbetas <- betas
            bestgammas <- gammas
            bestsigmas <- adm1.sigma
        }

        ## Perform a regression on each state
        stage1.betas <- matrix(NA, M, K)
        for (jj in 1:M) {
            included <- adm1 == jj
            if (sum(included) < K + 1)
                next

            dmxxsjj <- as.matrix(dmxxs[included,])
            modjj <- nnnpls(dmxxsjj, dmyy[included], stacked$coeff)

            ## Multiply back in exponent, if dmbins were generated with an assumed gamma
            for (kk in 1:K)
                stage1.betas[jj, kk] <- modjj$x[kk] * mean(dmxxsjj[, kk] / dmxxs.orig[included, kk], na.rm=T)
        }

        ## Prepare values for second stage regressions
        stage1.logbetas <- log(abs(stage1.betas))
        stage1.logbetas[stage1.logbetas == -Inf] <- NA

        ## Second stage regression for each predictor
        betas <- rep(NA, K)
        gammas <- rep(NA, sum(kls))

        gammas.so.far <- 0 # Keep track of how many coefficients used
        for (kk in 1:K) {
            valid <- is.finite(stage1.logbetas[, kk])

            gammas.here <- sum(kls[kk, ])
            if (gammas.here == 0)
                next

            modkk <- fastLm(cbind(rep(1, sum(valid)), zzs[valid, kls[kk, ]]), stage1.logbetas[valid, kk])

            ## Prepare values for log likelihood
            betas[kk] <- exp(modkk$coeff[1]) * sign(stacked$coeff[kk])
            gammas[(gammas.so.far+1):(gammas.so.far+gammas.here)] <- modkk$coeff[-1]
            gammas.so.far <- gammas.so.far + gammas.here
        }

        likeli2 <- calc.likeli.demeaned(dmxxs.orig, dmyy, zzs, kls, adm1, betas, gammas, adm1.sigma, weights)

        ## Report progress
        print(c(iter, log2(1/armijo.factor), as.integer(max(likeli1, likeli2))))

        ## Check if we have converged
        if (abs(likeli2 - bestlikeli) < 1e-6 || armijo.factor < 1e-6)
            break

        if (likeli2 > bestlikeli) {
            bestlikeli <- likeli2
            bestbetas <- betas
            bestgammas <- gammas
            bestsigmas <- adm1.sigma
        }

        if (max(likeli1, likeli2) < bestlikeli) {
            ## Oops!  We spun too far out
            gammas <- (bestgammas + armijo.factor * gammas) / (1 + armijo.factor)
            bestgravity <- T # Can't do beta and sigma calc yet, so leave to next iteration
            armijo.factor <- armijo.factor / 2
            ## Note: We may never achive as good a likelihood, since it may be inconsistent with the sigmas from stage1
        } else
            bestgravity <- F

        ## Setup for next iteration
        dmxxs <- calc.covariated.predictors(dmxxs.orig, zzs, kls, adm1, gammas)
        for (jj in 1:M) {
            included <- adm1 == jj
            dmyy.weighted[included] <- dmyy[included] / adm1.sigma[jj]
            dmxxs.weighted[included,] <- dmxxs[included,] / adm1.sigma[jj]
        }

        dmyy.weighted <- dmyy.weighted * sqrt(weights)
        for (kk in 1:K)
            dmxxs.weighted[, kk] <- dmxxs.weighted[, kk] * sqrt(weights)
    }

    list(betas=betas, gammas=gammas, sigmas=adm1.sigma)
}

## Perform a weighted regression with partialed-out predictors dmxxs
stacked.regression <- function(K, L, gammas, dmyy, dmxxs.orig, zzs, kls, mm, weights) {
    dmxxs <- calc.covariated.predictors(dmxxs.orig, zzs, kls, mm, gammas)
    for (kk in 1:K)
        dmxxs[, kk] <- dmxxs[, kk] * sqrt(weights)

    fastLm(dmxxs, dmyy * sqrt(weights))
}

## Get stacked betas
stacked.betas <- function(K, L, gammas, dmyy, dmxxs.orig, zzs, kls, mm, weights) {
    stacked <- stacked.regression(K, L, gammas, dmyy, dmxxs.orig, zzs, kls, mm, weights)

    stacked$coeff
}

## Estimate the maximum likelihood, using R's optim function
estimate.logspec.optim <- function(yy, xxs, zzs, kls, adm1, adm2, weights=1, initgammas=NULL) {
    list2env(check.arguments(yy, xxs, zzs, kls, adm1, adm2), environment())
    list2env(demean.yxs(K, yy, xxs, adm2), environment())
    dmxxs.orig <- dmxxs

    if (is.null(initgammas)) {
        ## Approximation 1: No covariate effect
        for (kk in 1:K)
            dmxxs[, kk] <- dmxxs[, kk] * sqrt(weights)
        stacked <- fastLm(dmxxs, dmyy * sqrt(weights))
        sigma <- sd(stacked$residuals)

        betas <- stacked$coeff
        gammas <- rep(0, sum(kls))
    } else {
        gammas <- initgammas
        stacked <- stacked.regression(K, L, gammas, dmyy, dmxxs.orig, zzs, kls, adm1, weights)
        betas <- stacked$coeff
        sigma <- sd(stacked$residuals)
    }

    print(c("Step 1:", calc.likeli.demeaned(dmxxs.orig, dmyy, zzs, kls, adm1, betas, gammas, rep(sigma, M), weights)))

    ## Approximation 2: Identically distributed
    objective <- function(params) {
        gammas <- params[1:sum(kls)]
        sigma <- rep(params[sum(kls)+1], M)

        betas <- stacked.betas(K, L, gammas, dmyy, dmxxs.orig, zzs, kls, adm1, weights)

        -calc.likeli.demeaned(dmxxs.orig, dmyy, zzs, kls, adm1, betas, gammas, sigma, weights)
    }

    params <- c(gammas, sigma)
    soln <- optim(params, objective)

    gammas <- soln$par[1:sum(kls)]
    sigma <- soln$par[sum(kls)+1]

    betas <- stacked.betas(K, L, gammas, dmyy, dmxxs.orig, zzs, kls, adm1, weights)

    print(c("Step 2:", calc.likeli.demeaned(dmxxs.orig, dmyy, zzs, kls, adm1, betas, gammas, rep(sigma, M), weights)))

    ## Approximation 3: State-clustered errors
    dmxxs <- calc.covariated.predictors(dmxxs.orig, zzs, kls, adm1, gammas)
    for (kk in 1:K)
        dmxxs[, kk] <- dmxxs[, kk] * sqrt(weights)
    stacked <- fastLm(dmxxs, dmyy * sqrt(weights))

    betas <- stacked$coeff

    sigma <- c()
    for (jj in 1:M)
        sigma <- c(sigma, sd(stacked$residuals[adm1 == jj]))

    print(c("Step 3:", calc.likeli.demeaned(dmxxs.orig, dmyy, zzs, kls, adm1, betas, gammas, sigma, weights)))

    objective2 <- function(params) {
        gammas <- params[1:sum(kls)]
        sigma <- params[(sum(kls)+1):length(params)]

        betas <- stacked.betas(K, L, gammas, dmyy, dmxxs.orig, zzs, kls, adm1, weights)

        -calc.likeli.demeaned(dmxxs.orig, dmyy, zzs, kls, adm1, betas, gammas, sigma, weights)
    }

    params <- c(soln$par[1:sum(kls)], sigma)
    soln <- optim(params, objective2)

    gammas <- soln$par[1:sum(kls)]
    sigma <- soln$par[(sum(kls)+1):length(soln$par)]

    betas <- stacked.betas(K, L, gammas, dmyy, dmxxs.orig, zzs, kls, adm1, weights)

    print(c("Step 4:", calc.likeli.demeaned(dmxxs.orig, dmyy, zzs, kls, adm1, betas, gammas, sigma, weights)))

    list(betas=betas, gammas=gammas, sigma=sigma)
}

source("methast.R")
source("ranges.R")
