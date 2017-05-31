source("logspec.R")

search.logspec <- function(yy, xxs, zzs, kls, adm1, adm2, weights=1, maxiter=100, initgammas=NULL) {
    list2env(check.arguments(yy, xxs, zzs, kls, adm1, adm2), environment())

    print("Preparing for search...")
    list2env(demean.yxs(K, yy, xxs, adm2), environment())

    search.logspec.demeaned(dmyy, dmxxs, zzs, kls, adm1, adm2, weights=weights, maxiter=maxiter, gammas=initgammas)
}

search.logspec.demeaned <- function(dmyy, dmxxs, zzs, kls, adm1, adm2, weights=1, maxiter=100, betas=NULL, gammas=NULL, sigmas=NULL, betases=NULL, gammases=NULL, bestlikeli=-Inf, skipmethod=0, eps=.01) {
    if (maxiter == 0)
        return(list(betas=betas, gammas=gammas, sigmas=sigmas, betases=betases, gammases=gammases, likeli=bestlikeli))

    list2env(check.arguments(dmyy, dmxxs, zzs, kls, adm1, adm2), environment())

    if (skipmethod != 1) {
        print("Search: Convergence")
        result <- estimate.logspec.demeaned(dmyy, dmxxs, zzs, kls, adm1, adm2, weights=weights, initgammas=gammas)
        if (result$likeli > bestlikeli) {
            print(paste("Improved!", result$likeli))

            if (result$likeli > bestlikeli + eps)
                return(search.logspec.demeaned(dmyy, dmxxs, zzs, kls, adm1, adm2, weights=weights, maxiter=maxiter-1, betas=result$betas, gammas=result$gammas, sigmas=result$sigmas, betases=betases, gammases=gammases, bestlikeli=result$likeli, skipmethod=1, eps=eps))
            else {
                betas <- result$betas
                gammas <- result$gammas
                sigmas <- result$sigmas
                bestlikeli <- result$likeli
            }
        }
    }

    if (skipmethod != 2) {
        print("Search: Gradient ascent")
        result <- estimate.logspec.gammaoptim.demeaned(dmyy, dmxxs, zzs, kls, adm1, adm2, sigmas=sigmas, weights=weights, initgammas=gammas)
        if (-result$value > bestlikeli) {
            print(paste("Improved!", -result$value))

            if (-result$value > bestlikeli + eps)
                return(search.logspec.demeaned(dmyy, dmxxs, zzs, kls, adm1, adm2, weights=weights, maxiter=maxiter-1, betas=NULL, gammas=result$par, sigmas=NULL, betases=betases, gammases=gammases, bestlikeli=-result$value, skipmethod=2, eps=eps))
            else {
                gammas <- result$par
                bestlikeli <- -result$value
            }
        }
    }

    if (skipmethod != 3) {
        print("Search: Optimization")
        result <- estimate.logspec.optim.demeaned(dmyy, dmxxs, zzs, kls, adm1, adm2, weights=weights, initgammas=gammas)
        likeli <- calc.likeli.demeaned(dmxxs, dmyy, zzs, kls, adm1, result$betas, result$gammas, result$sigma, weights)

        if (likeli > bestlikeli) {
            print(paste("Improved!", likeli))

            if (likeli > bestlikeli + eps)
                return(search.logspec.demeaned(dmyy, dmxxs, zzs, kls, adm1, adm2, weights=weights, maxiter=maxiter-1, betas=result$betas, gammas=result$gammas, sigmas=result$sigma, betases=betases, gammases=gammases, bestlikeli=likeli, skipmethod=3, eps=eps))
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
        if (is.null(betases)) {
            print("Estimating SEs")
            vcv <- calc.vcv.ols(K, L, dmxxs, dmyy, zzs, kls, adm1, betas, gammas, sigmas, weights)
            ses <- sqrt(abs(diag(vcv)))
            betases <- ses[1:K]
            gammases <- ses[(K+1):(K+sum(kls))]
        }

        result <- methast.betagamma(K, L, dmxxs, dmyy, zzs, kls, adm1, 600, betas, gammas, betases, gammases, weights=weights)
        if (result$best.likeli > bestlikeli) {
            print(paste("Improved!", result$best.likeli))

            ## Since only 1 seed, average with existing ses
            betases <- (betases + apply(result$betas[101:600,], 2, sd)) / 2
            gammases <- (gammases + apply(result$gammas[101:600,], 2, sd)) / 2

            print(betases)
            print(gammases)

            if (result$best.likeli > bestlikeli + eps)
                return(search.logspec.demeaned(dmyy, dmxxs, zzs, kls, adm1, adm2, weights=weights, maxiter=maxiter-1, betas=result$betas[result$best.index,], gammas=result$gammas[result$best.index,], sigmas=sigmas, betases=betases, gammases=gammases, bestlikeli=result$best.likeli, skipmethod=4, eps=eps))
            else {
                betas <- result$betas[result$best.index,]
                gammas <- result$gammas[result$best.index,]
                betases <- betases
                gammases <- gammases
                bestlikeli <- result$best.likeli
            }
        }
    }

    list(betas=betas, gammas=gammas, sigmas=sigmas, betases=betases, gammases=gammases, likeli=bestlikeli)
}
