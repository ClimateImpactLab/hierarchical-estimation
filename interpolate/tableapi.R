source("logspec.R")

ta.arguments <- function(df, outname, adm1name, adm2name, prednames, covarnames) {
    yy <- df[, outname]
    adm1 <- as.numeric(factor(df[, adm1name]))
    adm2 <- as.numeric(factor(paste(df[, adm1name], df[, adm2name], sep=', ')))

    unipreds <- unique(prednames)
    unicovars <- unique(covarnames)
    xxs <- df[, unipreds]

    ## Generate an order across representative rows
    zzs.representatives <- !duplicated(adm1)
    zzs.adm1 <- df[zzs.representatives, adm1name]
    zzs.adm1order <- c()
    for (mm in 1:max(adm1))
        zzs.adm1order <- c(zzs.adm1order, which(zzs.adm1 == mm))
    zzs <- df[zzs.representatives, unicovars][zzs.adm1order, ] # pull out, in order of adm1 numbers

    kls <- matrix(F, ncol(xxs), ncol(zzs))
    for (kk in 1:ncol(xxs))
        kls[kk, unicovars %in% covarnames[prednames == unipreds[kk]]] <- T

    list(yy=yy, adm1=adm1, adm2=adm2, xxs=xxs, zzs=zzs, kls=kls)
}

ta.estimate.logspec <- function(df, outname, adm1name, adm2name, prednames, covarnames, weights=1) {
    list2env(ta.arguments(df, outname, adm1name, adm2name, prednames, covarnames), environment())
    result <- estimate.logspec(yy, xxs, zzs, kls, adm1, adm2, weights=weights)
    print(result)
    estimate.logspec.optim(yy, xxs, zzs, kls, adm1, adm2, weights=weights, initgammas=result$gammas)
}

ta.estimate.vcv <- function(betas, gammas, sigmas, df, outname, adm1name, adm2name, prednames, covarnames, ...) {
    list2env(ta.arguments(df, outname, adm1name, adm2name, prednames, covarnames), environment())
    estimate.vcv(betas, gammas, sigmas, yy, xxs, zzs, kls, adm1, adm2, ...)
}
