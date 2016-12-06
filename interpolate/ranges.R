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

## Single Metropolis-Hastings with automatic tuning
## A 44% acceptance rate is optimal for 1-D sampling
## Want (a^44) * (b^56) = 1; Say a = 1.1, (1 / (1.1^44))^(1 / 56) = 0.9278487
methast <- function(K, L, dmxxs, dmyy, zzs, mm, iter, beta0, gamma0, betaerr, gammaerr) {
    params = matrix(NA, iter, K + K*L)
    params[1, ] = c(beta0, gamma0)

    last.likeli <- calc.likeli.nosigma(K, L, dmxxs, dmyy, zzs, mm, beta0, gamma0)

    ## Look out for an even better solution
    best.likeli <- last.likeli
    best.index <- 1

    sd.product <- 1

    for (ii in 2:iter) {
        if (ii %% 100 == 0)
            print(ii)
        beta.sample <- rnorm(K, params[ii-1, 1:K], betaerr * sd.product)
        gamma.sample <- matrix(rnorm(K*L, params[ii-1, (K+1):(K+K*L)], c(gammaerr) * sd.product), K, L)

        this.likeli <- calc.likeli.nosigma(K, L, dmxxs, dmyy, zzs, mm, beta.sample, gamma.sample)
        prob <- exp(this.likeli - last.likeli)

        if (min(prob, 1) > runif(1)) {
            params[ii, ] <- c(beta.sample, gamma.sample)
            last.likeli <- this.likeli
            sd.product <- sd.product * 1.1 # was too modest

            if (this.likeli > best.likeli) {
                best.likeli <- this.likeli
                best.index <- ii
            }
        } else {
            params[ii, ] <- params[ii-1, ]
            sd.product <- sd.product * 0.9278487 # was too bold
        }
    }

    list(betas=params[, 1:K], gammas=params[, (K+1):(K+K*L)], best.index=best.index)
}

## Use Metropolis-Hastings with N seeds
repeated.methast <- function(K, L, dmxxs, dmyy, zzs, adm1, iter, warmup, seeds, beta0, gamma0, betaerr, gammaerr) {
    betas <- matrix(NA, 0, K)
    gammas <- matrix(NA, 0, K*L)

    for (seed in 1:seeds) {
        print(paste("Seed", seed))
        methast.result <- methast(K, L, dmxxs, dmyy, zzs, adm1, iter, beta0, gamma0, betaerr, gammaerr)

        betas <- rbind(betas, methast.result$betas[(warmup+1):iter,])
        gammas <- rbind(gammas, methast.result$gammas[(warmup+1):iter,])

        if (methast.result$best.index != 1) {
            beta0 <- methast.result$betas[methast.result$best.index,]
            gamma0 <- as.matrix(methast.result$gammas[methast.result$best.index,], K, L)
        }
    }

    list(betas=betas, gammas=gammas, best.beta=beta0, best.gamma=gamma0)
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

estimate.vcv <- function(betas, gammas, sigmas, yy, xxs, zzs, adm1, adm2, iter=600, warmup=100, seed=4) {
    list2env(check.arguments(yy, xxs, zzs, adm1, adm2), parent.frame())

    ## De-mean observations
    dmyy <- regional.demean(yy, adm2)
    dmxxs <- xxs
    for (kk in 1:K)
        dmxxs[, kk] <- regional.demean(dmxxs[, kk], adm2)

    vcv.ols <- calc.vcv.ols(K, L, dmxxs, dmyy, zzs, adm1, betas, gammas, sigmas)
    se.ols <- sqrt(diag(vcv.ols))
    result <- repeated.methast(K, L, dmxxs, dmyy, zzs, adm1, iter, warmup, seed, betas, gammas, se.ols[1:K], se.ols[(K+1):(K+K*L)])

    serr <- serr.conservative(vcv.ols, cbind(result$betas, result$gammas))
    if (sum(serr != se.ols) == 0)
        list(betas=result$best.beta, gammas=result$best.gamma, vcv=vcv.ols, se=se.ols)
    else
        list(betas=result$best.beta, gammas=result$best.gamma, vcv=t(serr) %*% cor(params) %*% t(t(serr)), se=serr)
}
