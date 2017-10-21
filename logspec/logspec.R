library(nnls)

source("shared.R")

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
calc.expected.demeaned <- function(dmxxs, zzs, kls, mm, betas, gammas) {
    covariated.predictors <- calc.covariated.predictors(dmxxs, zzs, kls, mm, gammas)

    obsmean <- 0
    for (kk in 1:nrow(kls))
        obsmean <- obsmean + betas[kk] * covariated.predictors[, kk]

    obsmean
}

## Calculate the likelihood of the given parameters, after observations have been demeaned
calc.likeli.demeaned <- function(dmxxs, dmyy, zzs, kls, mm, betas, gammas, sigma, weights, prior) {
    obsmean <- calc.expected.demeaned(dmxxs, zzs, kls, mm, betas, gammas)
    sum(weights * dnorm(dmyy, obsmean, sigma[mm], log=T)) + prior(betas, gammas)
}

## Standard Gaussian gamma prior
gaussian.prior <- function(taus) {
    function(betas, gammas) {
        sum(dnorm(gammas, 0, taus, log=T))
    }
}

gaussian.gammapriorderiv <- function(taus) {
    function(gammas) {
        gammas / taus^2
    }
}

## No prior
noninformative.prior <- function(betas, gammas) {
    0
}

noninformative.gammapriorderiv <- function(gammas) {
    rep(0, length(gammas))
}

## Estimate the parameters of the log specification
estimate.logspec <- function(yy, xxs, zzs, kls, adm1, factors, weights=1, maxiter=1000, initgammas=NULL, prior=noninformative.prior, get.betas=stacked.betas) {
    list2env(check.arguments(yy, xxs, zzs, kls, adm1, factors), environment())
    list2env(demean.yxs(yy, xxs, factors, weights), environment())

    estimate.logspec.demeaned(dmyy, dmxxs, zzs, kls, adm1, weights=weights, maxiter=maxiter, initgammas=initgammas, prior=prior, get.betas=get.betas)
}

estimate.logspec.demeaned <- function(dmyy, dmxxs, zzs, kls, adm1, weights=1, maxiter=1000, initgammas=NULL, prior=noninformative.prior, get.betas=stacked.betas) {
    list2env(check.arguments(dmyy, dmxxs, zzs, kls, adm1), environment())

    if (maxiter > 1)
        print("Iterating...")

    bestlikeli <- -Inf # best likelihood we've seen
    bestgammas <- rep(0, sum(kls)) # gammas corresponding to bestlikeli
    bestbetas <- rep(0, K) # betas corresponding to bestlikeli
    bestsigmas <- rep(Inf, M) # adm1.sigma corresponding to bestlikeli
    armijo.factor <- 1 # the amount of movement away from bestgammas
    bestgravity <- F # should we move toward best this iteration

    if (is.null(initgammas))
        gammas <- bestgammas
    else
        gammas <- initgammas

    ## Will be weighted by region-error by dividing by standard error
    dmyy.weighted <- dmyy
    dmxxs.weighted <- dmxxs

    for (iter in 1:maxiter) {
        ## Perform stacked regression to get signs
        betas <- get.betas(K, L, gammas, dmyy.weighted, dmxxs.weighted, zzs, kls, adm1, weights)
        dmyy.weighted.hat <- calc.expected.demeaned(dmxxs.weighted, zzs, kls, adm1, betas, gammas)
        residuals <- dmyy.weighted - dmyy.weighted.hat

        ## Check the performance of this set of betas
        adm1.sigma <- rep(NA, M)
        for (jj in 1:M) {
            included <- adm1 == jj
            adm1.sigma[jj] <- mean(sd(residuals[included]) * dmyy[included] / dmyy.weighted[included])
        }
        adm1.sigma[is.na(adm1.sigma)] <- 1 # Dummy value (otherwise fails when all values in ADM1 identical)

        if (bestgravity) {
            betas <- (bestbetas + armijo.factor * betas) / (1 + armijo.factor)
            adm1.sigma <- (bestsigmas + armijo.factor * adm1.sigma) / (1 + armijo.factor)
        }

        likeli1 <- calc.likeli.demeaned(dmxxs, dmyy, zzs, kls, adm1, betas, gammas, adm1.sigma, weights, prior)

        ## Check if we have converged
        if (is.na(likeli1) || abs(likeli1 - bestlikeli) < 1e-6 || armijo.factor < 1e-6)
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
            modjj <- nnnpls(dmxxsjj, dmyy[included], betas)

            ## Multiply back in exponent, if dmbins were generated with an assumed gamma
            for (kk in 1:K)
                stage1.betas[jj, kk] <- modjj$x[kk]
        }

        ## Prepare values for second stage regressions
        stage1.logbetas <- log(abs(stage1.betas))
        stage1.logbetas[stage1.logbetas == -Inf] <- NA

        ## Second stage regression for each predictor-- only for gammas
        old.gammas <- gammas
        gammas <- rep(NA, sum(kls))

        gammas.so.far <- 0 # Keep track of how many coefficients used
        for (kk in 1:K) {
            valid <- is.finite(stage1.logbetas[, kk])

            gammas.here <- sum(kls[kk, ])
            if (gammas.here == 0)
                next

            if (sum(valid) < gammas.here) {
                gammas[(gammas.so.far+1):(gammas.so.far+gammas.here)] <- old.gammas[(gammas.so.far+1):(gammas.so.far+gammas.here)]
                gammas.so.far <- gammas.so.far + gammas.here
                next
            }

            yy.modkk <- stage1.logbetas[valid, kk]
            xx.modkk <- as.matrix(cbind(rep(1, sum(valid)), zzs[valid, kls[kk, ]]))
            modkk <- lm(yy.modkk ~ 0 + xx.modkk)

            ## Prepare values for log likelihood
            gammas[(gammas.so.far+1):(gammas.so.far+gammas.here)] <- modkk$coeff[-1]
            gammas.so.far <- gammas.so.far + gammas.here
        }

        ## Perform stacked regression to get betas
        betas <- get.betas(K, L, gammas, dmyy.weighted, dmxxs.weighted, zzs, kls, adm1, weights)

        likeli2 <- calc.likeli.demeaned(dmxxs, dmyy, zzs, kls, adm1, betas, gammas, adm1.sigma, weights, prior)

        ## Report progress
        print(paste(c(iter, log2(1/armijo.factor), format(max(likeli1, likeli2), digits=0, scientific=F)), collapse=" "))

        ## Check if we have converged
        if (is.na(likeli2) || abs(likeli2 - bestlikeli) < 1e-6 || armijo.factor < 1e-6)
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
        for (jj in 1:M) {
            included <- adm1 == jj
            dmyy.weighted[included] <- dmyy[included] / adm1.sigma[jj]
            dmxxs.weighted[included,] <- dmxxs[included,] / adm1.sigma[jj]
        }
    }

    betas <- get.betas(K, L, gammas, dmyy.weighted, dmxxs.weighted, zzs, kls, adm1, weights)

    list(betas=betas, gammas=gammas, sigmas=adm1.sigma, likeli=bestlikeli)
}

## Perform a weighted regression with partialed-out predictors dmxxs
stacked.regression <- function(K, L, gammas, dmyy, dmxxs.orig, zzs, kls, mm, weights) {
    dmxxs <- calc.covariated.predictors(dmxxs.orig, zzs, kls, mm, gammas)
    if (sum(!is.finite(dmxxs)) > 0)
        return(list(coeff=rep(NA, K)))

    if (length(weights) == 1)
        lm(dmyy ~ 0 + dmxxs)
    else
        lm(dmyy ~ 0 + dmxxs, weights=weights)
}

## Get betas that fit a stacked regression
stacked.betas <- function(K, L, gammas, dmyy, dmxxs.orig, zzs, kls, mm, weights) {
    stacked <- stacked.regression(K, L, gammas, dmyy, dmxxs.orig, zzs, kls, mm, weights)

    stacked$coeff
}

## Get betas that evaluate to a given betas0 at zz0
make.known.betas <- function(zz0, betas0) {
    function(K, L, gammas, dmyy, dmxxs.orig, zzs, kls, mm, weights) {
        expgamma0 <- rep(1, length(betas0)) # Only modify this for predictors with covariates
        gammas.so.far <- 0 # Keep track of how many coefficients used

        for (kk in 1:nrow(kls)) {
            gammas.here <- sum(kls[kk, ])
            if (gammas.here == 0)
                next # Nothing to do

            mygammas <- gammas[(gammas.so.far+1):(gammas.so.far+gammas.here)]
            gammas.so.far <- gammas.so.far + gammas.here
            expgamma0[kk] <- exp(sum(zz0[kls[kk, ]] * mygammas))
        }

        betas0 / expgamma0
    }
}

## Estimate the maximum likelihood, using R's optim function
estimate.logspec.optim <- function(yy, xxs, zzs, kls, adm1, factors, weights=1, initgammas=NULL, prior=noninformative.prior, get.betas=stacked.betas) {
    list2env(check.arguments(yy, xxs, zzs, kls, adm1, factors), environment())
    list2env(demean.yxs(yy, xxs, factors, weights), environment())

    estimate.logspec.optim.demeaned(dmyy, dmxxs, zzs, kls, adm1, weights=weights, initgammas=initgammas, prior=prior, get.betas=get.betas)
}

estimate.logspec.optim.demeaned <- function(dmyy, dmxxs, zzs, kls, adm1, weights=1, initgammas=NULL, prior=noninformative.prior, get.betas=stacked.betas) {
    list2env(check.arguments(dmyy, dmxxs, zzs, kls, adm1), environment())

    if (is.null(initgammas)) {
        ## Approximation 1: No covariate effect
        gammas <- rep(0, sum(kls))
    } else
        gammas <- initgammas

    betas <- get.betas(K, L, gammas, dmyy, dmxxs, zzs, kls, adm1, weights)

    dmyyhat <- calc.expected.demeaned(dmxxs, zzs, kls, adm1, betas, gammas)
    sigma <- sd(dmyy - dmyyhat)

    print(c("Step 1:", calc.likeli.demeaned(dmxxs, dmyy, zzs, kls, adm1, betas, gammas, rep(sigma, M), weights, prior)))

    ## Approximation 2: Identically distributed
    objective <- function(params) {
        gammas <- params[1:sum(kls)]
        sigma <- rep(params[sum(kls)+1], M)

        betas <- get.betas(K, L, gammas, dmyy, dmxxs, zzs, kls, adm1, weights)

        -calc.likeli.demeaned(dmxxs, dmyy, zzs, kls, adm1, betas, gammas, sigma, weights, prior)
    }

    params <- c(gammas, sigma)
    soln <- optim(params, objective)

    gammas <- soln$par[1:sum(kls)]
    sigma <- soln$par[sum(kls)+1]

    betas <- get.betas(K, L, gammas, dmyy, dmxxs, zzs, kls, adm1, weights)

    print(c("Step 2:", calc.likeli.demeaned(dmxxs, dmyy, zzs, kls, adm1, betas, gammas, rep(sigma, M), weights, prior)))

    ## Approximation 3: State-clustered errors
    dmyyhat <- calc.expected.demeaned(dmxxs, zzs, kls, adm1, betas, gammas)
    residuals <- dmyy - dmyyhat

    saved.sigma <- sigma

    sigma <- c()
    for (jj in 1:M)
        sigma <- c(sigma, sd(residuals[adm1 == jj]))

    ## Combine sigmas
    sigma <- (sigma + saved.sigma) / 2 # to handle sigmas of 0

    print(c("Step 3:", calc.likeli.demeaned(dmxxs, dmyy, zzs, kls, adm1, betas, gammas, sigma, weights, prior)))

    result <- estimate.logspec.optim.step4(dmyy, dmxxs, zzs, kls, adm1, soln$par[1:sum(kls)], sigma, weights=weights, prior=prior, get.betas=get.betas)

    print(c("Step 4:", calc.likeli.demeaned(dmxxs, dmyy, zzs, kls, adm1, betas, gammas, sigma, weights, prior)))

    result
}

estimate.logspec.optim.step4 <- function(dmyy, dmxxs, zzs, kls, adm1, gammas, sigma, weights=1, prior=noninformative.prior, get.betas=stacked.betas) {
    list2env(check.arguments(dmyy, dmxxs, zzs, kls, adm1), environment())

    objective2 <- function(params) {
        gammas <- params[1:sum(kls)]
        sigma <- params[(sum(kls)+1):length(params)]

        betas <- get.betas(K, L, gammas, dmyy, dmxxs, zzs, kls, adm1, weights)

        -calc.likeli.demeaned(dmxxs, dmyy, zzs, kls, adm1, betas, gammas, sigma, weights, prior)
    }

    params <- c(gammas, sigma)
    soln <- optim(params, objective2)

    gammas <- soln$par[1:sum(kls)]
    sigma <- soln$par[(sum(kls)+1):length(soln$par)]

    betas <- get.betas(K, L, gammas, dmyy, dmxxs, zzs, kls, adm1, weights)

    list(betas=betas, gammas=gammas, sigma=sigma)
}

source("methast.R")
source("ranges.R")
source("gradient.R")
source("predict.R")
