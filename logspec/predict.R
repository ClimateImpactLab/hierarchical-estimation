logspec.predict <- function(xxs, zzs, kls, adm1, betas, gammas, means=NULL) {
    list2env(check.arguments(rep(0, nrow(xxs)), xxs, zzs, kls, adm1), environment())

    yy <- calc.expected.demeaned(xxs, zzs, kls, adm1, betas, gammas)

    if (!is.null(means))
        yy <- yy + means

    yy
}

logspec.predict.betas <- function(zzs, kls, betas, gammas) {
    result <- t(matrix(betas, length(betas), nrow(zzs))) # Only modify this for predictors with covariates

    for (kk in 1:nrow(kls)) {
        if (any(kls[kk, ] > 0))
            next # Nothing to do: already beta

        result[, kk] <- betas[kk] * exp(as.matrix(zzs[, kls[kk, ], drop=F]) %*% gammas[kls[kk, ]])
    }

    result
}

rsqr.demeaned <- function(dmyy, dmxxs, zzs, kls, adm1, betas, gammas, weights=1) {
    dmyy.pred <- calc.expected.demeaned(dmxxs, zzs, kls, adm1, betas, gammas) ## NOTE: should include weights, when change by ADM1

    1 - sum(weights * (dmyy - dmyy.pred)^2) / sum(weights * dmyy^2)
}

rsqr <- function(yy, xxs, zzs, kls, adm1, factors, betas, gammas, weights=1) {
    means <- logspec.get.fe(yy, xxs, zzs, kls, adm1, factors, betas, gammas, weights=weights)
    yy.pred <- logspec.predict(xxs, zzs, kls, adm1, betas, gammas, means)

    if (length(weights) == 1)
        return(1 - sum((yy - yy.pred)^2) / sum((yy - mean(yy))^2))
    else
        return(1 - sum(weights * (yy - yy.pred)^2) / sum(weights * (yy - weighted.mean(yy, weights))^2))
}

rsqr.projected <- function(yy, xxs, zzs, kls, adm1, factors, betas, gammas, weights=1) {
    list2env(check.arguments(yy, xxs, zzs, kls, adm1, factors), environment())
    list2env(demean.yxs(yy, xxs, factors, weights), environment())

    rsqr.demeaned(dmyy, dmxxs, zzs, kls, adm1, betas, gammas, weights)
}
