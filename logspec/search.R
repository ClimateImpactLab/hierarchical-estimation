source("logspec.R")

search.logspec <- function(yy, xxs, zzs, kls, adm1, factors, weights=1, maxiter=100, initgammas=NULL, prior=noninformative.prior, gammapriorderiv=noninformative.gammapriorderiv, known.betas.info=NULL) {
    list2env(check.arguments(yy, xxs, zzs, kls, adm1, factors), environment())

    print("Preparing for search...")
    list2env(demean.yxs(yy, xxs, factors, weights), environment())

    search.logspec.demeaned(dmyy, dmxxs, zzs, kls, adm1, weights=weights, maxiter=maxiter, gammas=initgammas, prior=prior, gammapriorderiv=gammapriorderiv, known.betas.info=known.betas.info)
}

search.logspec.demeaned <- function(dmyy, dmxxs, zzs, kls, adm1, weights=1, maxiter=100, betas=NULL, gammas=NULL, sigmas=NULL, betases=NULL, gammases=NULL, bestlikeli=-Inf, skipmethod=0, eps=.01, prior=noninformative.prior, gammapriorderiv=noninformative.gammapriorderiv, known.betas.info=NULL) {
    if (maxiter == 0)
        return(list(betas=betas, gammas=gammas, sigmas=sigmas, betases=betases, gammases=gammases, likeli=bestlikeli))

    list2env(check.arguments(dmyy, dmxxs, zzs, kls, adm1), environment())

    if (is.null(known.betas.info))
        get.betas <- stacked.betas
    else
        get.betas <- make.known.betas(known.betas.info$zz0, known.betas.info$betas0)

    if (skipmethod != 1) {
        print("Search: Convergence")
        result <- estimate.logspec.demeaned(dmyy, dmxxs, zzs, kls, adm1, weights=weights, initgammas=gammas, prior=prior, get.betas=get.betas)
        if (result$likeli > bestlikeli) {
            print(paste("Improved!", result$likeli))

            if (result$likeli > bestlikeli + eps)
                return(search.logspec.demeaned(dmyy, dmxxs, zzs, kls, adm1, weights=weights, maxiter=maxiter-1, betas=result$betas, gammas=result$gammas, sigmas=result$sigmas, betases=betases, gammases=gammases, bestlikeli=result$likeli, skipmethod=1, eps=eps, prior=prior, known.betas.info=known.betas.info))
            else {
                betas <- result$betas
                gammas <- result$gammas
                sigmas <- result$sigmas
                bestlikeli <- result$likeli
            }
        }
    }

    if (skipmethod != 2 && sum(kls) > 0 && !is.null(sigmas)) {
        print("Search: Gradient ascent")
        result <- estimate.logspec.gammaoptim.demeaned(dmyy, dmxxs, zzs, kls, adm1, sigmas=sigmas, weights=weights, initgammas=gammas, prior=prior, gammapriorderiv=gammapriorderiv, get.betas=get.betas)
        if (-result$value > bestlikeli) {
            print(paste("Improved!", -result$value))

            betas <- get.betas(K, L, result$par, dmyy, dmxxs, zzs, kls, adm1, weights)

            if (-result$value > bestlikeli + eps)
                return(search.logspec.demeaned(dmyy, dmxxs, zzs, kls, adm1, weights=weights, maxiter=maxiter-1, betas=betas, gammas=result$par, sigmas=sigmas, betases=betases, gammases=gammases, bestlikeli=-result$value, skipmethod=2, eps=eps, prior=prior, known.betas.info=known.betas.info))
            else {
                gammas <- result$par
                bestlikeli <- -result$value
            }
        }
    }

    if (skipmethod != 3) {
        print("Search: Optimization")
        if (is.null(sigmas))
            result <- estimate.logspec.optim.demeaned(dmyy, dmxxs, zzs, kls, adm1, weights=weights, initgammas=gammas, prior=prior)
        else
            result <- estimate.logspec.optim.step4(dmyy, dmxxs, zzs, kls, adm1, gammas, sigmas, weights=weights, prior=prior, get.betas=get.betas)

        likeli <- calc.likeli.demeaned(dmxxs, dmyy, zzs, kls, adm1, result$betas, result$gammas, result$sigma, weights, prior)

        if (likeli > bestlikeli) {
            print(paste("Improved!", likeli))

            if (likeli > bestlikeli + eps)
                return(search.logspec.demeaned(dmyy, dmxxs, zzs, kls, adm1, weights=weights, maxiter=maxiter-1, betas=result$betas, gammas=result$gammas, sigmas=result$sigma, betases=betases, gammases=gammases, bestlikeli=likeli, skipmethod=3, eps=eps, prior=prior, known.betas.info=known.betas.info))
            else {
                betas <- result$betas
                gammas <- result$gammas
                sigmas <- result$sigma
                bestlikeli <- likeli
            }
        }
    }

    if (skipmethod != 4) {
        print("Search: Metropolis-Hastings")
        if (is.null(gammases)) {
            print("Estimating SEs")
            if (is.null(known.betas.info)) {
                ses <- tryCatch({
                    vcv <- calc.vcv.ols(K, L, dmxxs, dmyy, zzs, kls, adm1, betas, gammas, sigmas, weights, prior=prior)
                    sqrt(abs(diag(vcv)))
                }, error=function(e) {
                    result.each <- repeated.methast.each(K, L, dmxxs, dmyy, zzs, kls, adm1, 200, 100, 4, betas, gammas, sigmas, weights, prior=prior, verbose=F)
                    c(result.each$betaerr, result.each$gammaerr)
                })
                betases <- ses[1:K]
                gammases <- ses[(K+1):(K+sum(kls))]
            } else {
                gammases <- tryCatch({
                    vcv <- calc.vcv.ols.gammaonly(K, L, dmxxs, dmyy, zzs, kls, adm1, gammas, sigmas, weights, prior=prior)
                    sqrt(abs(diag(vcv)))
                }, error=function(e) {
                    repeated.methast.each.gammaonly(K, L, dmxxs, dmyy, zzs, kls, adm1, 200, 100, 4, betas, gammas, sigmas, weights, prior=prior, verbose=F)
                })
            }
        }

        print("Performing MCMC")
        if (is.null(known.betas.info)) {
            result <- methast.betagamma.sigma(K, L, dmxxs, dmyy, zzs, kls, adm1, 600, betas, gammas, sigmas, betases, gammases, weights=weights, prior=prior)

            betases.new <- apply(result$betas[101:600,], 2, sd)
        } else {
            result <- methast.gamma.sigma(K, L, dmxxs, dmyy, zzs, kls, adm1, 600, get.betas, gammas, sigmas, gammases, weights=weights, prior=prior)

            betases.new <- NA
        }

        if (result$best.likeli > bestlikeli) {
            print(paste("Improved!", result$best.likeli))

            ## Since only 1 seed, average with existing ses
            betases <- (betases + betases.new) / 2
            gammases <- (gammases + apply(result$gammas[101:600,], 2, sd)) / 2

            if (result$best.likeli > bestlikeli + eps)
                return(search.logspec.demeaned(dmyy, dmxxs, zzs, kls, adm1, weights=weights, maxiter=maxiter-1, betas=result$betas[result$best.index,], gammas=result$gammas[result$best.index,], sigmas=sigmas, betases=betases, gammases=gammases, bestlikeli=result$best.likeli, skipmethod=4, eps=eps, prior=prior, known.betas.info=known.betas.info))
            else {
                betas <- result$betas[result$best.index,]
                gammas <- result$gammas[result$best.index,]
                betases <- betases
                gammases <- gammases
                bestlikeli <- result$best.likeli
            }
        }
    }

    if (skipmethod != 5 && sum(kls) > 0) {
        print("Recentering covariates")

        ## Decide on offsets
        zzs.offset <- matrix(NA, nrow(zzs), 0)
        offsets <- c()
        for (ll in 1:L) {
            offset <- mean(zzs[, ll])
            if (offset / sd(zzs[, ll]) < 1e-3)
                offset <- rnorm(1, 0, sd(zzs[, ll])) # Choose offset at random

            zzs.offset <- cbind(zzs.offset, zzs[, ll] - offset)
            offsets <- c(offsets, offset)
        }

        if (is.null(known.betas.info))
            known.betas.info.offset <- NULL
        else
            known.betas.info.offset <- list(zz0=known.betas.info$zz0 - offsets, betas0=known.betas.info$betas0)

        result <- search.logspec.demeaned(dmyy, dmxxs, zzs.offset, kls, adm1, weights=weights, maxiter=1, betas=zoffset.adjust.betas(betas, offsets, kls, gammas), gammas=gammas, sigmas=sigmas, betases=zoffset.adjust.betas(betases, offsets, kls, gammas), gammases=gammases, bestlikeli=bestlikeli, skipmethod=5, eps=eps, prior=prior, known.betas.info=known.betas.info.offset)

        ## Re-calculate with true zzs
        if (is.null(known.betas.info))
            readjusted.betas <- zoffset.adjust.betas(result$betas, -offsets, kls, result$gammas)
        else
            readjusted.betas <- get.betas(K, L, result$gammas, dmyy, dmxxs, zzs, kls, adm1, weights)

        likeli <- calc.likeli.demeaned(dmxxs, dmyy, zzs, kls, adm1, readjusted.betas, result$gammas, result$sigma, weights, prior)

        if (likeli > bestlikeli) {
            print(paste("Improved!", likeli))

            readjusted.betases <- zoffset.adjust.betas(result$betases, -offsets, kls, result$gammas)
            if (likeli > bestlikeli + eps)
                return(search.logspec.demeaned(dmyy, dmxxs, zzs, kls, adm1, weights=weights, maxiter=maxiter-1, betas=readjusted.betas, gammas=result$gammas, sigmas=result$sigmas, betases=readjusted.betases, gammases=result$gammases, bestlikeli=likeli, skipmethod=5, eps=eps, prior=prior, known.betas.info=known.betas.info))
            else {
                betas <- readjusted.betas
                gammas <- result$gammas
                sigmas <- result$sigmas
                betases <- readjusted.betases
                gammases <- result$gammases
                bestlikeli <- likeli
            }
        }
    }

    list(betas=betas, gammas=gammas, sigmas=sigmas, betases=betases, gammases=gammases, likeli=bestlikeli)
}

zoffset.adjust.betas <- function(betas, offsets, kls, gammas) {
    gammas.so.far <- 0 # Keep track of how many coefficients used

    for (kk in 1:nrow(kls)) {
        gammas.here <- sum(kls[kk, ])
        if (gammas.here == 0)
            next # Nothing to do

        mygammas <- gammas[(gammas.so.far+1):(gammas.so.far+gammas.here)]
        gammas.so.far <- gammas.so.far + gammas.here
        betas[kk] <- betas[kk] * exp(sum(offsets[kls[kk, ]] * mygammas))
    }

    betas
}
