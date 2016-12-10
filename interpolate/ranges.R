## This library supports logspec.R, and logspec.R must be loaded.

## Calculate the log likelihood, computing ADM1 sigmas from residuals
calc.likeli.nosigma <- function(K, L, dmxxs, dmyy, zzs, mm, betas, gammas) {
    dmyy.exp <- calc.expected.demeaned(K, L, dmxxs, dmyy, zzs, mm, betas, gammas)

    sigmas <- c()
    for (jj in unique(mm)) {
        included <- mm == jj
        residuals <- dmyy.exp[included] - dmyy[included]
        sigmas <- c(sigmas, sd(residuals))
    }

    calc.likeli.demeaned(K, L, dmxxs, dmyy, zzs, mm, betas, gammas, sigmas)
}

make.methast.betagamma.likeli <- function(K, L, dmxxs, dmyy, zzs, mm) {
    function(param) {
        beta <- param[1:K]
        gamma <- param[(K+1):(K*K*L)]
        calc.likeli.nosigma(K, L, dmxxs, dmyy, zzs, mm, beta, matrix(gamma, K, L))
    }
}

methast.betagamma <- function(K, L, dmxxs, dmyy, zzs, mm, iter, beta0, gamma0, betaerr, gammaerr) {
    result <- methast(iter, c(beta0, gamma0), c(betaerr, gammaerr),
                      make.methast.betagamma.likeli(K, L, dmxxs, dmyy, zzs, mm))

    list(betas=result$params[, 1:K], gammas=result$params[, (K+1):(K+K*L)], best.index=result$best.index)
}

## Use Metropolis-Hastings with N seeds
repeated.methast.betagamma <- function(K, L, dmxxs, dmyy, zzs, adm1, iter, warmup, seeds, beta0, gamma0, betaerr, gammaerr) {
    result <- repeated.methast(seeds, iter, warmup,
                               c(beta0, gamma0), c(betaerr, gammaerr),
                               make.methast.betagamma.likeli(K, L, dmxxs, dmyy, zzs, adm1))

    list(betas=result$params[, 1:K], gammas=result$params[, (K+1):(K+K*L)], best.beta=result$best.param[1:K], best.gamma=as.matrix(result$best.param[(K+1):(K+K*L)], K, L))
}

## Single Metropolis-Hastings with automatic tuning
repeated.methast.each <- function(K, L, dmxxs, dmyy, zzs, mm, iter, warmup, seeds, beta0, gamma0, sigmas) {
    betaerr <- c()
    for (kk in 1:length(beta0)) {
        result <- repeated.methast(seeds, iter, warmup, beta0[kk], 1,
                                   function(beta) {
                                       beta2 = beta0
                                       beta2[kk] <- beta
                                       calc.likeli.demeaned(K, L, dmxxs, dmyy, zzs, mm,
                                                            beta2, gamma0, sigmas)
                                   })
        betaerr <- c(betaerr, sd(result$params))
    }

    gammaerr <- c()
    for (kk in 1:dim(gamma0)[1]) {
        for (ll in 1:dim(gamma0)[2]) {
            result <- repeated.methast(seeds, iter, warmup, gamma0[kk, ll], 1,
                                       function(gamma) {
                                           gamma2 = gamma0
                                           gamma2[kk, ll] <- gamma
                                           calc.likeli.demeaned(K, L, dmxxs, dmyy, zzs, mm,
                                                                beta0, gamma2, sigmas)
                                       })
            gammaerr <- c(gammaerr, sd(result$params))

        }
    }
    gammaerr <- matrix(gammaerr, K, L)

    list(betaerr=betaerr, gammaerr=gammaerr)
}

calc.vcv.ols <- function(K, L, dmxxs, dmyy, zzs, adm1, betas, gammas, sigmas) {
    objective <- function(params) {
        betas <- params[1:K]
        gammas <- matrix(params[(K+1):((L+1)*K)], K, L)
        -calc.likeli.demeaned(K, L, dmxxs, dmyy, zzs, adm1, betas, gammas, sigmas)
    }

    params <- c(betas, as.vector(gammas))
    soln <- optim(params, objective, hessian=T)

    solve(soln$hessian)
}

calc.vcv.methast <- function(K, L, dmxxs, dmyy, zzs, adm1, iter, warmup, seeds, betas, gammas, sigmas) {
    result <- repeated.methast.each(K, L, dmxxs, dmyy, zzs, adm1, iter, warmup, seeds, betas, gammas, sigmas)

    vcv.bayes(K, L, dmxxs, dmyy, zzs, adm1, iter, warmup, seeds, betas, gammas, result$betaerr, result$gammaerr)
}

vcv.bayes <- function(K, L, dmxxs, dmyy, zzs, adm1, iters, warmup, seeds, beta0, gamma0, betaerr, gammaerr) {
    result <- repeated.methast(K, L, dmxxs, dmyy, zzs, adm1, iter, warmup, seeds, beta0, gamm0, betaerr, gammaerr)

    cov(cbind(result$betas, result$gammas))
}

ses.tails <- function(params) {
    ## Also, calculate effective SD that would fit the 95% range
    obsd <- apply(apply(params, 2, function(xx) quantile(xx, probs=c(.025, .975))), 2, diff)
    expd <- qnorm(.975) - qnorm(.025)
    obsd / expd
}

serr.conservative <- function(vcv.ols, params) {
    sd.ols <- sqrt(diag(vcv.ols))
    sd.bayes <- apply(params, 2, sd)
    sd.tails <- ses.tails(params)

    pmax(sd.ols, sd.bayes, sd.tails)
}

estimate.vcv <- function(betas, gammas, sigmas, yy, xxs, zzs, adm1, adm2, iter=600, warmup=100, seeds=4, use.ols=T) {
    list2env(check.arguments(yy, xxs, zzs, adm1, adm2), parent.frame())
    list2env(demean.yxs(yy, xxs, adm2), parent.frame())

    if (use.ols) {
        vcv.start <- tryCatch({
            calc.vcv.ols(K, L, dmxxs, dmyy, zzs, adm1, betas, gammas, sigmas)
        }, error=function(e) {
            NULL
        })

        if (is.null(vcv.start))
            use.ols <- F
    }

    if (use.ols) {
        se.start <- sqrt(diag(vcv.start))
    } else {
        result.each <- repeated.methast.each(K, L, dmxxs, dmyy, zzs, adm1, iter, warmup, seeds, betas, gammas, sigmas)
        se.start <- c(result.each$betaerr, result.each$gammaerr)
    }

    result <- repeated.methast.betagamma(K, L, dmxxs, dmyy, zzs, adm1, iter, warmup, seeds, betas, gammas, se.start[1:K], matrix(se.start[(K+1):(K+K*L)], K, L))
    if (!use.ols)
        vcv.start <- cov(cbind(result$betas, result$gammas))

    serr <- serr.conservative(vcv.start, cbind(result$betas, result$gammas))
    if (sum(serr != se.start) == 0 && use.ols)
        list(betas=result$best.beta, gammas=result$best.gamma, vcv=vcv.start, se=se.start)
    else
        list(betas=result$best.beta, gammas=result$best.gamma, vcv=t(serr) %*% cor(cbind(result$betas, result$gammas)) %*% t(t(serr)), se=serr)
}

estimate.se <- function(betas, gammas, sigmas, yy, xxs, zzs, adm1, adm2, iter=600, warmup=100, seeds=4, use.ols=T) {
    list2env(check.arguments(yy, xxs, zzs, adm1, adm2), parent.frame())
    list2env(demean.yxs(yy, xxs, adm2), parent.frame())

    if (use.ols) {
        vcv.start <- tryCatch({
            calc.vcv.ols(K, L, dmxxs, dmyy, zzs, adm1, betas, gammas, sigmas)
        }, error=function(e) {
            NULL
        })

        if (is.null(vcv.start))
            use.ols <- F
    }

    if (use.ols) {
        return(sqrt(diag(vcv.start)))
    } else {
        result.each <- repeated.methast.each(K, L, dmxxs, dmyy, zzs, adm1, iter, warmup, seeds, betas, gammas, sigmas)
        return(c(result.each$betaerr, result.each$gammaerr))
    }
}
