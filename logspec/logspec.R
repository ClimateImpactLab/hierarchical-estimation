library(nnls)

source("shared.R")

calc.covariated.predictors <- function(dmxxs, zzs, kls, mm, gammas) {
    ##' Calculate x_k exp(sum gamma_kl z_l) as a NxK matrix
    ##' @param dmxxs NxK matrix of predictors
    ##' @param zzs MxL matrix of covariates
    ##' @param kls KxL matrix of 0 - G; 0 means exclude covariate l on predictor k, >0 means use gamma G
    ##' @param mm N vector of 1 - M for the region
    ##' @param gammas G vector for the set of gammas

    result <- dmxxs # Only modify this for predictors with covariates

    for (kk in 1:nrow(kls)) {
        if (all(kls[kk, ] == 0))
            next # Nothing to do: already dmxxs[, kk]

        result[, kk] <- dmxxs[, kk] * exp(as.matrix(zzs[, kls[kk, ] > 0]) %*% gammas[kls[kk, ]])[mm]
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
    bestgammas <- rep(0, max(kls)) # gammas corresponding to bestlikeli
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

        ## Second stage regression for each predictor-- only for gammas
        old.gammas <- gammas
        gammas <- rep(NA, max(kls))

        if (max(kls) == L * sum(rowSums(kls) > 0)) { # different gammas for each covariate
            ## Prepare values for second stage regressions
            stage1.logbetas <- log(abs(stage1.betas))
            stage1.logbetas[stage1.logbetas == -Inf] <- NA

            for (kk in 1:K) {
                valid <- is.finite(stage1.logbetas[, kk])

                if (all(kls[kk, ] == 0))
                    next

                if (sum(valid) < sum(kls[kk, ] > 0)) {
                    gammas[kls[kk, ]] <- old.gammas[kls[kk, ]]
                    next
                }

                yy.modkk <- stage1.logbetas[valid, kk]
                xx.modkk <- as.matrix(cbind(rep(1, sum(valid)), zzs[valid, kls[kk, ] > 0]))
                modkk <- lm(yy.modkk ~ 0 + xx.modkk)

                ## Prepare values for log likelihood
                gammas[kls[kk, ]] <- modkk$coeff[-1]
            }
        } else {
            uniquekls <- unique(kls)
            klssigs <- apply(kls, 1, function(row) paste(row, collapse=","))

            for (rr in 1:nrow(uniquekls)) {
                ## Find all k that share this uniquekls
                uniqueklssig <- paste(uniquekls[rr,], collapse=",")
                kincluded <- klssigs == uniqueklssig

                stage1.logbx <- log(abs(rowSums(stage1.betas[adm1, kincluded] * dmxxs[, kincluded])))
                valid <- is.finite(stage1.logbx)

                yy.mod2 <- stage1.logbx[valid]
                XX.mod2 <- as.matrix(rep(1, sum(valid)), sum(valid), 1)
                for (gg in uniquekls[rr, uniquekls[rr,] > 0])
                    XX.mod2 <- cbind(XX.mod2, zzs[adm1[valid], uniquekls[rr,] == gg])

                modkk <- lm(yy.mod2 ~ 0 + XX.mod2)

                ## Prepare values for log likelihood
                gammas[uniquekls[rr, ]] <- modkk$coeff[-1]
            }
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

        for (kk in 1:nrow(kls)) {
            if (all(kls[kk, ] == 0))
                next # Nothing to do

            expgamma0[kk] <- exp(sum(zz0[kls[kk, ] > 0] * gammas[kls[kk, ]]))
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
        gammas <- rep(0, max(kls))
    } else
        gammas <- initgammas

    betas <- get.betas(K, L, gammas, dmyy, dmxxs, zzs, kls, adm1, weights)

    dmyyhat <- calc.expected.demeaned(dmxxs, zzs, kls, adm1, betas, gammas)
    sigma <- sd(dmyy - dmyyhat)

    print(c("Step 1:", calc.likeli.demeaned(dmxxs, dmyy, zzs, kls, adm1, betas, gammas, rep(sigma, M), weights, prior)))

    ## Approximation 2: Identically distributed
    objective <- function(params) {
        gammas <- params[1:max(kls)]
        sigma <- rep(params[max(kls)+1], M)

        betas <- get.betas(K, L, gammas, dmyy, dmxxs, zzs, kls, adm1, weights)

        -calc.likeli.demeaned(dmxxs, dmyy, zzs, kls, adm1, betas, gammas, sigma, weights, prior)
    }

    params <- c(gammas, sigma)
    soln <- optim(params, objective)

    gammas <- soln$par[1:max(kls)]
    sigma <- soln$par[max(kls)+1]

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

    result <- estimate.logspec.optim.step4(dmyy, dmxxs, zzs, kls, adm1, soln$par[1:max(kls)], sigma, weights=weights, prior=prior, get.betas=get.betas)

    print(c("Step 4:", calc.likeli.demeaned(dmxxs, dmyy, zzs, kls, adm1, betas, gammas, sigma, weights, prior)))

    result
}

estimate.logspec.optim.step1.knownsigma <- function(dmyy, dmxxs, zzs, kls, adm1, gammas, sigma, weights=1, prior=noninformative.prior, get.betas=stacked.betas) {
    list2env(check.arguments(dmyy, dmxxs, zzs, kls, adm1), environment())

    objective1 <- function(gammas) {
        betas <- get.betas(K, L, gammas, dmyy, dmxxs, zzs, kls, adm1, weights)
        -calc.likeli.demeaned(dmxxs, dmyy, zzs, kls, adm1, betas, gammas, sigma, weights, prior)
    }

    soln <- optim(gammas, objective1)

    betas <- get.betas(K, L, soln$par, dmyy, dmxxs, zzs, kls, adm1, weights)

    list(betas=betas, gammas=soln$par)
}

estimate.logspec.optim.step4 <- function(dmyy, dmxxs, zzs, kls, adm1, gammas, sigma, weights=1, prior=noninformative.prior, get.betas=stacked.betas) {
    list2env(check.arguments(dmyy, dmxxs, zzs, kls, adm1), environment())

    objective2 <- function(params) {
        gammas <- params[1:max(kls)]
        sigma <- params[(max(kls)+1):length(params)]

        betas <- get.betas(K, L, gammas, dmyy, dmxxs, zzs, kls, adm1, weights)

        -calc.likeli.demeaned(dmxxs, dmyy, zzs, kls, adm1, betas, gammas, sigma, weights, prior)
    }

    params <- c(gammas, sigma)
    soln <- optim(params, objective2)

    gammas <- soln$par[1:max(kls)]
    sigma <- soln$par[(max(kls)+1):length(soln$par)]

    betas <- get.betas(K, L, gammas, dmyy, dmxxs, zzs, kls, adm1, weights)

    list(betas=betas, gammas=gammas, sigma=sigma)
}

source("methast.R")
source("ranges.R")
source("gradient.R")
source("predict.R")
