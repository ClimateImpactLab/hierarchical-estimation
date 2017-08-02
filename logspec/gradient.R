calc.gamma.gradient <- function(dmxxs, dmyy, zzs, kls, mm, betas, gammas, sigmas, weights) {
    obsmean <- calc.expected.demeaned(dmxxs, zzs, kls, mm, betas, gammas)

    klcombos <- which(kls, arr.ind=T)

    term1.all <- 2 * (dmyy - obsmean)
    term2.byk <- calc.covariated.predictors(dmxxs, zzs, kls, mm, gammas)
    term3.bym <- (-1 / (2 * sigmas^2))[mm]

    gradients <- rep(NA, nrow(klcombos))

    for (ii in 1:nrow(klcombos)) {
        kk <- klcombos[ii, 1]
        ll <- klcombos[ii, 2]

        gradients[ii] <- sum(term1.all * betas[kk] * term2.byk[, kk] * zzs[mm, ll] * weights * term3.bym)
    }

    return(gradients)
}

get.gamma.gradient <- function(yy, xxs, zzs, kls, adm1, adm2, betas, gammas, sigmas, weights=1) {
    list2env(check.arguments(yy, xxs, zzs, kls, adm1, adm2), environment())
    list2env(demean.yxs(K, yy, xxs, adm2, weights), environment())

    calc.gamma.gradient(dmxxs, dmyy, zzs, kls, adm1, betas, gammas, sigmas, weights)
}

## This uses numerical differentiation; it's a good comparison
estimate.logspec.gammaoptim.nograd <- function(yy, xxs, zzs, kls, adm1, adm2, sigmas, weights=1, initgammas=NULL, prior=noninformative.prior, get.betas=stacked.betas) {
    list2env(check.arguments(yy, xxs, zzs, kls, adm1, adm2), environment())
    list2env(demean.yxs(K, yy, xxs, adm2, weights), environment())

    if (is.null(initgammas))
        initgammas <- rep(0, sum(kls))

    objective <- function(gammas) {
        betas <- get.betas(K, L, gammas, dmyy, dmxxs, zzs, kls, adm1, weights)

        -calc.likeli.demeaned(dmxxs, dmyy, zzs, kls, adm1, betas, gammas, sigmas, weights, prior)
    }

    optim(initgammas, objective)
}

# The same as above, but with gradients
estimate.logspec.gammaoptim <- function(yy, xxs, zzs, kls, adm1, adm2, sigmas, weights=1, initgammas=NULL, prior=noninformative.prior, gammapriorderiv=noninformative.gammapriorderiv, get.betas=stacked.betas) {
    list2env(check.arguments(yy, xxs, zzs, kls, adm1, adm2), environment())
    list2env(demean.yxs(K, yy, xxs, adm2, weights), environment())

    estimate.logspec.gammaoptim.demeaned(dmyy, dmxxs, zzs, kls, adm1, adm2, sigmas, weights=weights, initgammas=initgammas, prior=prior, gammapriorderiv=gammapriorderiv, get.betas=get.betas)
}

estimate.logspec.gammaoptim.demeaned <- function(dmyy, dmxxs, zzs, kls, adm1, adm2, sigmas, weights=1, initgammas=NULL, prior=noninformative.prior, gammapriorderiv=noninformative.gammapriorderiv, get.betas=stacked.betas) {
    list2env(check.arguments(dmyy, dmxxs, zzs, kls, adm1, adm2), environment())
    if (is.null(initgammas))
        initgammas <- rep(0, sum(kls))

    objective <- function(gammas) {
        betas <- get.betas(K, L, gammas, dmyy, dmxxs, zzs, kls, adm1, weights)

        -calc.likeli.demeaned(dmxxs, dmyy, zzs, kls, adm1, betas, gammas, sigmas, weights, prior)
    }

    gradient <- function(gammas) {
        betas <- get.betas(K, L, gammas, dmyy, dmxxs, zzs, kls, adm1, weights)

        calc.gamma.gradient(dmxxs, dmyy, zzs, kls, adm1, betas, gammas, sigmas, weights) + gammapriorderiv(gammas)
    }

    optim(initgammas, objective, gradient, method="BFGS")
}

## Determine \partial^2 f / \partial x \partial z, the marginal effect
## of the interaction between x and z at given covariates
marginal.interactions <- function(zz, kls, betas, gammas) {
    marginals <- c()

    for (kk in 1:length(betas)) {
        gammas.here <- sum(kls[kk, ])
        if (gammas.here == 0)
            next # Nothing to do: already dmxxs[, kk]

        gammas.so.far <- length(marginals)
        mygammas <- gammas[(gammas.so.far+1):(gammas.so.far+gammas.here)]

        mymarginals <- betas[kk] * mygammas * exp(sum(zz[kls[kk, ]] * mygammas))
        marginals <- c(marginals, mymarginals)
    }

    marginals
}
