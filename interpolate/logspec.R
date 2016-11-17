library(nnls)
library(RcppArmadillo)

calc.expected.demeaned <- function(K, dmxxs, dmyy, zzs, mm, betas, gammas) {
    obsmean <- 0
    for (kk in 1:K)
        obsmean <- obsmean + betas[kk] * dmxxs[, kk] * exp(as.matrix(zzs[, ((kk-1)*L + 1):(kk*L)]) %*% gammas[kk, ])[mm]

    obsmean
}

calc.likeli.demeaned <- function(K, dmxxs, dmyy, zzs, mm, betas, gammas, sigma) {
    obsmean <- calc.expected.demeaned(K, dmxxs, dmyy, zzs, mm, betas, gammas)
    sum(dnorm(dmyy, obsmean, sigma[mm], log=T))
}

regional.demean <- function(values, regions) {
    for (region in unique(regions)) {
        regioniis <- which(regions == region)
        values[regioniis] <- values[regioniis] - mean(values[regioniis])
    }

    values
}

estimate.logspec <- function(yy, xxs, zzs, adm1, adm2, maxiter=1000) {
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

    dmyy <- regional.demean(yy, adm2)

    ## Demean all predictors
    dmxxs.orig <- xxs
    for (kk in 1:K)
        dmxxs.orig[, kk] <- regional.demean(xxs[, kk], adm2)

    ## Start with no covariate effects
    dmxxs <- dmxxs.orig

    print("Iterating...")

    bestlikeli <- -Inf # best likelihood we've seen
    bestgammas <- rep(0, K*L) # gammas corresponding to bestlikeli
    armijo.factor <- 1 # the amount of movement away from bestgammas

    for (iter in 1:maxiter) {
        ## Perform stacked regression to get signs
        stacked <- fastLm(dmxxs, dmyy)

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
            modkk <- fastLm(cbind(rep(1, M), zzs[, ((kk-1)*L + 1):(kk*L)]), stage1.logbetas[, kk])

            ## Prepare values for log likelihood
            betas[kk] <- exp(modkk$coeff[1]) * sign(stacked$coeff[kk])
            gammas[kk, ] <- modkk$coeff[-1]
        }

        likeli <- calc.likeli.demeaned(K, dmxxs.orig, dmyy, zzs, adm1, betas, gammas, stage1.sigma)
        ## Report progress
        print(c(iter, likeli))

        ## Check if we have converged
        if (abs(likeli - bestlikeli) < 1e-6 || armijo.factor < 1e-6)
            break

        if (likeli > bestlikeli) {
            bestlikeli <- likeli
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
    }

    print("Calculating Hessian...")

    ## Estimate the Hessian
    objective <- function(params) {
        betas <- params[1:K]
        gammas <- matrix(params[(K+1):((L+1)*K)], K, L)
        -calc.likeli.demeaned(K, dmxxs.orig, dmyy, zzs, adm1, betas, gammas, stage1.sigma)
    }

    params <- c(betas, as.vector(gammas))
    soln <- optim(params, objective, hessian=T)

    ses <- sqrt(abs(diag(solve(soln$hessian))))

    list(betas=soln$par[1:K], gammas=soln$par[(K+1):((L+1)*K)], ses.betas=ses[1:K], ses.gammas=ses[(K+1):((L+1)*K)])
}
