## This library supports logspec.R, and logspec.R must be loaded.

## Calculate the log likelihood, computing ADM1 sigmas from residuals
calc.likeli.nosigma <- function(dmxxs, dmyy, zzs, kls, mm, betas, gammas, weights, prior) {
    dmyy.exp <- calc.expected.demeaned(dmxxs, zzs, kls, mm, betas, gammas)

    sigmas <- c()
    for (jj in unique(mm)) {
        included <- mm == jj
        residuals <- dmyy.exp[included] - dmyy[included]
        sigmas <- c(sigmas, sd(residuals))
    }

    calc.likeli.demeaned(dmxxs, dmyy, zzs, kls, mm, betas, gammas, sigmas, weights, prior)
}

calc.likeli.withsigma <- function(dmxxs, dmyy, zzs, kls, mm, betas, gammas, sigmas, weights, prior) {
    calc.likeli.demeaned(dmxxs, dmyy, zzs, kls, mm, betas, gammas, sigmas, weights, prior)
}

make.methast.betagamma.likeli <- function(K, L, dmxxs, dmyy, zzs, kls, mm, weights, prior) {
    function(param) {
        beta <- param[1:K]
        gamma <- param[(K+1):(K+sum(kls))]
        calc.likeli.nosigma(dmxxs, dmyy, zzs, kls, mm, beta, gamma, weights, prior)
    }
}

make.methast.betagamma.sigma.likeli <- function(K, L, dmxxs, dmyy, zzs, kls, mm, sigmas, weights, prior) {
    function(param) {
        betas <- param[1:K]
        gammas <- param[(K+1):(K+sum(kls))]
        calc.likeli.demeaned(dmxxs, dmyy, zzs, kls, mm, betas, gammas, sigmas, weights, prior)
    }
}

make.methast.gamma.sigma.likeli <- function(K, L, dmxxs, dmyy, zzs, kls, mm, sigmas, weights, prior, get.betas) {
    function(gammas) {
        betas <- get.betas(K, L, gammas, dmyy, dmxxs.orig, zzs, kls, mm, weights)
        calc.likeli.demeaned(dmxxs, dmyy, zzs, kls, mm, betas, gammas, sigmas, weights, prior)
    }
}

methast.betagamma <- function(K, L, dmxxs, dmyy, zzs, kls, mm, iter, beta0, gamma0, betaerr, gammaerr, weights=1, prior=noninformative.prior) {
    likelifunc <- make.methast.betagamma.likeli(K, L, dmxxs, dmyy, zzs, kls, mm, weights, prior)
    result <- methast(iter, c(beta0, gamma0), c(betaerr, gammaerr), likelifunc)

    list(betas=result$params[, 1:K], gammas=result$params[, (K+1):(K+sum(kls))], best.likeli=likelifunc(result$params[result$best.index,]), best.index=result$best.index)
}

methast.betagamma.sigma <- function(K, L, dmxxs, dmyy, zzs, kls, mm, iter, beta0, gamma0, sigmas, betaerr, gammaerr, weights=1, prior=noninformative.prior) {
    likelifunc <- make.methast.betagamma.sigma.likeli(K, L, dmxxs, dmyy, zzs, kls, mm, sigmas, weights, prior)
    result <- methast(iter, c(beta0, gamma0), c(betaerr, gammaerr), likelifunc)

    if (sum(kls) > 0)
        list(betas=result$params[, 1:K], gammas=result$params[, (K+1):(K+sum(kls))], best.likeli=likelifunc(result$params[result$best.index,]), best.index=result$best.index)
    else
        list(betas=result$params[, 1:K], gammas=matrix(NA, nrow(result$params), 0), best.likeli=likelifunc(result$params[result$best.index,]), best.index=result$best.index)
}

methast.gamma.sigma <- function(K, L, dmxxs, dmyy, zzs, kls, mm, iter, get.betas, gamma0, sigmas, gammaerr, weights=1, prior=noninformative.prior) {
    likelifunc <- make.methast.gamma.sigma.likeli(K, L, dmxxs, dmyy, zzs, kls, mm, sigmas, weights, prior, get.betas)

    if (sum(kls) == 0) { # Nothing to do
        betas <- get.betas(K, L, gammas, dmyy, dmxxs, zzs, kls, mm, weights)
        return(list(betas=betas, gammas=gamma0, best.likeli=likelifunc(gamma0), best.index=NA))
    }

    result <- methast(iter, gamma0, gammaerr, likelifunc)
    betas0 <- get.betas(K, L, result$params[result$best.index, ], dmyy, dmxxs, zzs, kls, mm, weights)
    betas <- t(matrix(betas0, length(betas0), iter))

    list(betas=betas, gammas=result$params, best.likeli=likelifunc(result$params[result$best.index,]), best.index=result$best.index)
}

## Use Metropolis-Hastings with N seeds
repeated.methast.betagamma <- function(K, L, dmxxs, dmyy, zzs, kls, adm1, iter, warmup, seeds, beta0, gamma0, betaerr, gammaerr, weights=1, prior=noninformative.prior) {
    result <- repeated.methast(seeds, iter, warmup,
                               c(beta0, gamma0), c(betaerr, gammaerr),
                               make.methast.betagamma.likeli(K, L, dmxxs, dmyy, zzs, kls, adm1, weights, prior))

    list(betas=result$params[, 1:K], gammas=result$params[, (K+1):(K+sum(kls))], best.beta=result$best.param[1:K], best.gamma=result$best.param[(K+1):(K+sum(kls))])
}

## Single Metropolis-Hastings with automatic tuning
repeated.methast.each <- function(K, L, dmxxs, dmyy, zzs, kls, mm, iter, warmup, seeds, beta0, gamma0, sigmas, weights=1, prior=noninformative.prior, verbose=F) {
    betaerr <- c()
    for (kk in 1:length(beta0)) {
        print(c("Beta", kk))
        result <- repeated.methast(seeds, iter, warmup, beta0[kk], 1,
                                   function(beta) {
                                       beta2 = beta0
                                       beta2[kk] <- beta
                                       calc.likeli.demeaned(dmxxs, dmyy, zzs, kls, mm,
                                                            beta2, gamma0, sigmas, weights, prior)
                                   }, verbose=verbose)
        betaerr <- c(betaerr, sd(result$params))
    }

    gammaerr <- repeated.methast.each.gammaonly(K, L, dmxxs, dmyy, zzs, kls, mm, iter, warmup, seeds, beta0, gamma0, sigmas, weights, prior, verbose)

    list(betaerr=betaerr, gammaerr=gammaerr)
}

repeated.methast.each.gammaonly <- function(K, L, dmxxs, dmyy, zzs, kls, mm, iter, warmup, seeds, beta0, gamma0, sigmas, weights=1, prior=noninformative.prior, verbose=F) {
    gammaerr <- c()
    for (kl in 1:length(gamma0)) {
        print(c("Gamma", kl))
        result <- repeated.methast(seeds, iter, warmup, gamma0[kl], 1,
                                   function(gamma) {
                                       gamma2 = gamma0
                                       gamma2[kl] <- gamma
                                       calc.likeli.demeaned(dmxxs, dmyy, zzs, kls, mm,
                                                            beta0, gamma2, sigmas, weights, prior)
                                   }, verbose=verbose)
        gammaerr <- c(gammaerr, sd(result$params))
    }

    gammaerr
}

calc.vcv.ols <- function(K, L, dmxxs, dmyy, zzs, kls, adm1, betas, gammas, sigmas, weights, prior=noninformative.prior) {
    objective <- function(params) {
        betas <- params[1:K]
        gammas <- params[(K+1):(K+sum(kls))]
        -calc.likeli.demeaned(dmxxs, dmyy, zzs, kls, adm1, betas, gammas, sigmas, weights, prior)
    }

    params <- c(betas, gammas)
    soln <- optim(params, objective, hessian=T)

    solve(soln$hessian)
}

calc.vcv.ols.gammaonly <- function(K, L, dmxxs, dmyy, zzs, kls, adm1, gammas, sigmas, weights, prior=noninformative.prior, get.betas=stacked.betas) {
    objective <- function(gammas) {
        betas <- get.betas(K, L, gammas, dmyy, dmxxs, zzs, kls, adm1, weights)
        -calc.likeli.demeaned(dmxxs, dmyy, zzs, kls, adm1, betas, gammas, sigmas, weights, prior)
    }

    soln <- optim(gammas, objective, hessian=T)

    solve(soln$hessian)
}

calc.vcv.methast <- function(K, L, dmxxs, dmyy, zzs, kls, adm1, iter, warmup, seeds, betas, gammas, sigmas, weights=1, prior=noninformative.prior) {
    result <- repeated.methast.each(K, L, dmxxs, dmyy, zzs, kls, adm1, iter, warmup, seeds, betas, gammas, sigmas, weights=1, prior=prior)

    vcv.bayes(K, L, dmxxs, dmyy, zzs, kls, adm1, iter, warmup, seeds, betas, gammas, result$betaerr, result$gammaerr, weights=1)
}

vcv.bayes <- function(K, L, dmxxs, dmyy, zzs, kls, adm1, iters, warmup, seeds, beta0, gamma0, betaerr, gammaerr, weights=1) {
    result <- repeated.methast(K, L, dmxxs, dmyy, zzs, kls, adm1, iter, warmup, seeds, beta0, gamma0, betaerr, gammaerr, weights=1)

    cov(cbind(result$betas, result$gammas))
}

ses.tails <- function(params) {
    ## Also, calculate effective SD that would fit the 95% range
    obsd <- apply(apply(params, 2, function(xx) quantile(xx, probs=c(.025, .975))), 2, diff)
    expd <- qnorm(.975) - qnorm(.025)
    obsd / expd
}

serr.conservative <- function(vcv.ols, params) {
    sd.ols <- sqrt(abs(diag(vcv.ols))) # Can get small negative values
    sd.bayes <- apply(params, 2, sd)
    sd.tails <- ses.tails(params)

    pmax(sd.ols, sd.bayes, sd.tails)
}

estimate.vcv <- function(betas, gammas, sigmas, yy, xxs, zzs, kls, adm1, adm2, iter=600, warmup=100, seeds=4, use.ols=T, weights=1, prior=noninformative.prior) {
    list2env(check.arguments(yy, xxs, zzs, kls, adm1, adm2), environment())
    list2env(demean.yxs(K, yy, xxs, adm2, weights), environment())

    if (use.ols) {
        vcv.start <- tryCatch({
            calc.vcv.ols(K, L, dmxxs, dmyy, zzs, kls, adm1, betas, gammas, sigmas, weights, prior=prior)
        }, error=function(e) {
            NULL
        })

        if (is.null(vcv.start))
            use.ols <- F
    }

    if (use.ols) {
        se.start <- sqrt(abs(diag(vcv.start))) # Can get small negative values
    } else {
        result.each <- repeated.methast.each(K, L, dmxxs, dmyy, zzs, kls, adm1, iter, warmup, seeds, betas, gammas, sigmas, weights=1, prior=prior)
        se.start <- c(result.each$betaerr, result.each$gammaerr)
    }

    result <- repeated.methast.betagamma(K, L, dmxxs, dmyy, zzs, kls, adm1, iter, warmup, seeds, betas, gammas, se.start[1:K], se.start[(K+1):(K+sum(kls))], weights=1, prior=prior)
    if (!use.ols)
        vcv.start <- cov(cbind(result$betas, result$gammas))

    serr <- serr.conservative(vcv.start, cbind(result$betas, result$gammas))
    if (sum(serr != se.start) == 0 && use.ols)
        list(betas=result$best.beta, gammas=result$best.gamma, vcv=vcv.start, se=se.start)
    else
        list(betas=result$best.beta, gammas=result$best.gamma, vcv=diag(serr) %*% cor(cbind(result$betas, result$gammas)) %*% diag(serr), se=serr)
}

estimate.se <- function(betas, gammas, sigmas, yy, xxs, zzs, kls, adm1, adm2, iter=600, warmup=100, seeds=4, use.ols=T, weights=1, prior=noninformative.prior) {
    list2env(check.arguments(yy, xxs, zzs, kls, adm1, adm2), environment())
    list2env(demean.yxs(K, yy, xxs, adm2, weights), environment())

    if (use.ols) {
        vcv.start <- tryCatch({
            calc.vcv.ols(K, L, dmxxs, dmyy, zzs, kls, adm1, betas, gammas, sigmas, weights=1, prior=prior)
        }, error=function(e) {
            NULL
        })

        if (is.null(vcv.start))
            use.ols <- F
    }

    if (use.ols) {
        return(sqrt(abs(diag(vcv.start)))) # Can get small negative values
    } else {
        result.each <- repeated.methast.each(K, L, dmxxs, dmyy, zzs, kls, adm1, iter, warmup, seeds, betas, gammas, sigmas, weights=1, prior=prior)
        return(c(result.each$betaerr, result.each$gammaerr))
    }
}
