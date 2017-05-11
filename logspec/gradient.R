calc.gamma.gradient <- function(dmxxs, dmyy, zzs, kls, mm, betas, gammas, sigmas, weights) {
    obsmean <- calc.expected.demeaned(dmxxs, dmyy, zzs, kls, mm, betas, gammas)

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
    list2env(demean.yxs(K, yy, xxs, adm2), environment())

    calc.gamma.gradient(dmxxs, dmyy, zzs, kls, adm1, betas, gammas, sigmas, weights)
}

## This uses numerical differentiation; it's a good comparison
estimate.logspec.gammaoptim.nograd <- function(yy, xxs, zzs, kls, adm1, adm2, sigmas, weights=1, initgammas=NULL) {
    list2env(check.arguments(yy, xxs, zzs, kls, adm1, adm2), environment())
    list2env(demean.yxs(K, yy, xxs, adm2), environment())

    if (is.null(initgammas))
        initgammas <- rep(0, sum(kls))

    objective <- function(gammas) {
        betas <- stacked.betas(K, L, gammas, dmyy, dmxxs, zzs, kls, adm1, weights)

        -calc.likeli.demeaned(dmxxs, dmyy, zzs, kls, adm1, betas, gammas, sigmas, weights)
    }

    optim(initgammas, objective)
}

# The same as above, but with gradients
estimate.logspec.gammaoptim <- function(yy, xxs, zzs, kls, adm1, adm2, sigmas, weights=1, initgammas=NULL) {
    list2env(check.arguments(yy, xxs, zzs, kls, adm1, adm2), environment())
    list2env(demean.yxs(K, yy, xxs, adm2), environment())

    estimate.logspec.gammaoptim.demeaned(dmyy, dmxxs, zzs, kls, adm1, adm2, sigmas, weights=weights, initgammas=initgammas)
}

estimate.logspec.gammaoptim.demeaned <- function(dmyy, dmxxs, zzs, kls, adm1, adm2, sigmas, weights=1, initgammas=NULL) {
    list2env(check.arguments(dmyy, dmxxs, zzs, kls, adm1, adm2), environment())
    if (is.null(initgammas))
        initgammas <- rep(0, sum(kls))

    objective <- function(gammas) {
        betas <- stacked.betas(K, L, gammas, dmyy, dmxxs, zzs, kls, adm1, weights)

        -calc.likeli.demeaned(dmxxs, dmyy, zzs, kls, adm1, betas, gammas, sigmas, weights)
    }

    gradient <- function(gammas) {
        betas <- stacked.betas(K, L, gammas, dmyy, dmxxs, zzs, kls, adm1, weights)

        gradient <- calc.gamma.gradient(dmxxs, dmyy, zzs, kls, adm1, betas, gammas, sigmas, weights)
        gradient
    }

    optim(initgammas, objective, gradient, method="BFGS")
}