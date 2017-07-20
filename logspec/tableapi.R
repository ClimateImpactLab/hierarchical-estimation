source("logspec.R")
source("search.R")

ta.arguments <- function(df, outname, adm1name, adm2name, prednames, covarnames) {
    yy <- df[, outname]

    adm1 <- as.numeric(factor(df[, adm1name]))
    adm2 <- as.numeric(factor(paste(df[, adm1name], df[, adm2name], sep=', ')))

    unipreds <- unique(prednames)
    unicovars <- unique(covarnames)
    if (!('1' %in% unicovars))
        stop("One of the covariates must be '1' for each predictor.")
    unicovars <- unicovars[unicovars != '1'] # Drop this covariate

    xxs <- df[, unipreds]

    ## Generate an order across representative rows
    zzs.representatives <- !duplicated(adm1)
    zzs.adm1 <- df[zzs.representatives, adm1name]
    zzs.adm1order <- c()
    for (mm in 1:max(adm1))
        zzs.adm1order <- c(zzs.adm1order, which(zzs.adm1 == mm))
    zzs <- df[zzs.representatives, unicovars, drop=F][zzs.adm1order, , drop=F] # pull out, in order of adm1 numbers

    kls <- matrix(F, ncol(xxs), ncol(zzs))
    for (kk in 1:ncol(xxs))
        kls[kk, unicovars %in% covarnames[prednames == unipreds[kk]]] <- T

    list(yy=yy, adm1=adm1, adm2=adm2, xxs=xxs, zzs=zzs, kls=kls)
}

## Wrapper on estimate.logspec
ta.estimate.logspec <- function(df, outname, adm1name, adm2name, prednames, covarnames, weights=1) {
    list2env(ta.arguments(df, outname, adm1name, adm2name, prednames, covarnames), environment())

    search.logspec(yy, xxs, zzs, kls, adm1, adm2, weights=weights)
}

## Wrapper on estimate.vcv
ta.estimate.vcv <- function(betas, gammas, sigmas, df, outname, adm1name, adm2name, prednames, covarnames, ...) {
    list2env(ta.arguments(df, outname, adm1name, adm2name, prednames, covarnames), environment())
    estimate.vcv(betas, gammas, sigmas, yy, xxs, zzs, kls, adm1, adm2, ...)
}

ta.rsqr <- function(fit, df, outname, adm1name, adm2name, prednames, covarnames) {
    list2env(ta.arguments(df, outname, adm1name, adm2name, prednames, covarnames), environment())

    rsqr(yy, xxs, zzs, kls, adm1, adm2, fit$betas, fit$gammas)
}
