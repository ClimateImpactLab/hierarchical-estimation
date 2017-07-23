logspec.get.fe <- function(yy, xxs, zzs, kls, adm1, adm2, betas, gammas) {
    list2env(check.arguments(yy, xxs, zzs, kls, adm1, adm2), environment())

    pred.yy <- calc.expected.demeaned(xxs, zzs, kls, adm1, betas, gammas) # Okay that not demeaned here

    fes <- list()
    for (region in unique(adm2)) {
        regioniis <- which(adm2 == region)
        fes[[region]] <- mean(yy[regioniis] - pred.yy[regioniis])
    }

    fes
}

logspec.predict <- function(xxs, zzs, kls, adm1, adm2, betas, gammas, fes=NULL) {
    list2env(check.arguments(rep(0, nrow(xxs)), xxs, zzs, kls, adm1, adm2), environment())

    yy <- calc.expected.demeaned(xxs, zzs, kls, adm1, betas, gammas)

    if (!is.null(fes)) {
        for (region in unique(adm2)) {
            regioniis <- which(adm2 == region)
            yy[regioniis] <- yy[regioniis] + fes[[region]]
        }
    }

    yy
}

rsqr.demeaned <- function(dmyy, dmxxs, zzs, kls, adm1, adm2, betas, gammas, weights=1) {
    dmyy.pred <- calc.expected.demeaned(dmxxs, zzs, kls, adm1, betas, gammas) ## NOTE: should include weights, when change by ADM1

    1 - sum(weights * (dmyy - dmyy.pred)^2) / sum(weights * dmyy^2)
}

rsqr <- function(yy, xxs, zzs, kls, adm1, adm2, betas, gammas, weights=1) {
    fes <- logspec.get.fe(yy, xxs, zzs, kls, adm1, adm2, betas, gammas)
    yy.pred <- logspec.predict(xxs, zzs, kls, adm1, adm2, betas, gammas, fes)

    if (length(weights) == 1)
        return(1 - sum((yy - yy.pred)^2) / sum((yy - mean(yy))^2))
    else
        return(1 - sum(weights * (yy - yy.pred)^2) / sum(weights * (yy - weighted.mean(yy, weights))^2))
}

rsqr.projected <- function(yy, xxs, zzs, kls, adm1, adm2, betas, gammas, weights=1) {
    list2env(check.arguments(yy, xxs, zzs, kls, adm1, adm2), environment())
    list2env(demean.yxs(K, yy, xxs, adm2), environment())

    rsqr.demeaned(dmyy, dmxxs, zzs, kls, adm1, adm2, betas, gammas, weights)
}
