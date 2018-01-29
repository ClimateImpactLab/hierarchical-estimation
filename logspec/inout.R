## Follows https://docs.google.com/document/d/1V7Wc4r-sgREwl70bAHv19tUrA6DS-w8IAPUMTS3Sy_8/edit
ta.output.csvv <- function(df, outname, adm1name, prednames, covarnames, factorouts, fit, vcvres, oneline, version, dependencies, description, vardefs, varunits, filepath, weights=1, demeanfile=NULL) {
    require(yaml)
    list2env(ta.arguments(df, outname, adm1name, prednames, covarnames, factorouts, demeanfile), environment())
    list2env(check.arguments(yy, xxs, zzs, kls, adm1, factors), environment())
    list2env(demean.yxs(yy, xxs, factors, weights), environment())

    fp <- file(filepath, 'w')

    variables <- list(outcome=paste0(vardefs[[outname]], " [", varunits[[outname]], "]"))
    for (name in unique(c(prednames, covarnames))) {
        if (name == '1')
            next
        variables[[name]] <- paste0(vardefs[[name]], " [", varunits[[name]], "]")
    }

    cat("---\n", file=fp)
    cat(as.yaml(list(oneline=oneline, version=version, dependencies=paste(dependencies, collapse=", "), description=description,
                     "csvv-version"="girdin-2017-01-10", variables=variables)), file=fp)
    cat("...\n", file=fp)

    cat("observations\n", file=fp)
    cat(nrow(df), file=fp)
    cat("\nprednames\n", file=fp)
    cat(paste(prednames, collapse=','), file=fp)
    cat("\ncovarnames\n", file=fp)
    cat(paste(covarnames, collapse=','), file=fp)

    cat("\ngamma\n", file=fp)

    orderedgammas <- c()
    vcvorder <- c()
    betas.so.far <- 0
    gammas.so.far <- 0
    for (ii in 1:length(prednames)) {
        if (covarnames[ii] == '1') {
            betas.so.far <- betas.so.far + 1
            orderedgammas <- c(orderedgammas, fit$betas[betas.so.far])
            vcvorder <- c(vcvorder, betas.so.far)
        } else {
            gammas.so.far <- gammas.so.far + 1
            orderedgammas <- c(orderedgammas, fit$gammas[gammas.so.far])
            vcvorder <- c(vcvorder, length(fit$betas) + gammas.so.far)
        }
    }

    cat(paste(orderedgammas, collapse=','), file=fp)

    cat("\ngammavcv\n", file=fp)
    for (ii in length(vcvorder))
        cat(paste(vcvfit$vcv[vcvorder[ii], vcvorder], collapse=','), file=fp)

    ## Predict, so that can calculate residual variance
    dmyyhat <- logspec.predict(dmxxs, zzs, kls, adm1, fit$betas, fit$gammas)
    cat("\nresidvcv\n", file=fp)
    cat(sum((dmyy - dmyyhat)^2) / (nrow(df) - 2), file=fp)
    cat("\n", file=fp)
    close(fp)
}
