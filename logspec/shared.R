library(lfe)

## Check that all of the given arguments are correctly specified
check.arguments <- function(yy, xxs, zzs, kls, adm1, factors=NULL) {
    if (is.null(nrow(xxs)) || is.null(ncol(xxs)))
        stop("xxs must be a matrix or data.frame.")

    if (is.null(nrow(zzs)) || is.null(ncol(zzs)))
        stop("zzs must be a matrix or data.frame.")

    N <- length(yy)
    if (nrow(xxs) != N || length(adm1) != N)
        stop("yy, xxs, and adm1 must all have the same number of observations.")

    if (!is.null(factors)) {
        if (is.list(factors)) {
            for (ii in 1:length(factors))
                if (length(factors[[ii]]) != N)
                    stop("All factors must the same number of observations as yy.")
        } else
            if (length(factors) != N)
                stop("ADM2 must the same number of observations as yy.")
    }

    K <- ncol(xxs)
    L <- ncol(zzs)

    if (nrow(kls) != K)
        stop("kls must have as many rows as xxs has columns.")

    if (ncol(kls) != L)
        stop("kls must have as many columns as zzs.")

    if (!is.logical(kls))
        stop("kls must consist only of trues and falses.")

    M <- max(adm1)
    if (length(unique(adm1)) != M)
        stop("adm1 must contain intgers from 1 to M.")

    if (nrow(zzs) != M)
        stop("zzs must have the same number of rows as ADM1 values.")

    return(list(K=K, L=L, N=N, M=M))
}

## Demean a set of values by region (partial out region fixed effects)
regional.demean <- function(values, regions, weights) {
    if (length(weights) == 1) {
        for (region in unique(regions)) {
            regioniis <- which(regions == region)
            values[regioniis] <- values[regioniis] - mean(values[regioniis])
        }
    } else {
        for (region in unique(regions)) {
            regioniis <- which(regions == region)
            values[regioniis] <- values[regioniis] - weighted.mean(values[regioniis], weights[regioniis])
        }
    }

    values
}

demean.yxs <- function(yy, xxs, factors, weights) {
    if (is.list(factors))
        demean.yxs.lfe(yy, xxs, factors, weights)
    else
        demean.yxs.adm(yy, xxs, factors, weights)
}

## Demean all observations by ADM2 regions (partial out ADM2 fixed effects)
demean.yxs.adm <- function(yy, xxs, adm2, weights) {
    ## De-mean observations
    dmyy <- regional.demean(yy, adm2, weights)
    dmxxs <- xxs
    for (kk in 1:ncol(dmxxs))
        dmxxs[, kk] <- regional.demean(dmxxs[, kk], adm2, weights)

    list(dmyy=dmyy, dmxxs=as.matrix(dmxxs))
}

## Demean, based on a factor list
demean.yxs.lfe <- function(yy, xxs, fl, weights) {
    results <- demeanlist(list(yy, xxs), fl, weights=weights)

    list(dmyy=results[[1]], dmxxs=as.matrix(results[[2]]))
}

logspec.get.fe <- function(yy, xxs, zzs, kls, adm1, factors, betas, gammas, weights=1) {
    list2env(check.arguments(yy, xxs, zzs, kls, adm1, factors), environment())

    pred.yy <- calc.expected.demeaned(xxs, zzs, kls, adm1, betas, gammas) # Okay that not demeaned here

    means <- rep(NA, length(yy))

    if (is.list(factors)) {
        means <- demeanlist(list(yy, xxs), factors, weights=weights, means=T)
    } else {
        for (region in unique(factors)) {
            regioniis <- which(factors == region)
            if (length(weights) == 1)
                means[regioniis] <- mean(yy[regioniis] - pred.yy[regioniis])
            else
                means[regioniis] <- weighted.mean(yy[regioniis] - pred.yy[regioniis], weights[regioniis])
        }
    }

    means
}
