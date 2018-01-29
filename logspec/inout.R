## Follows https://docs.google.com/document/d/1V7Wc4r-sgREwl70bAHv19tUrA6DS-w8IAPUMTS3Sy_8/edit
ta.output.csvv <- function(df, outname, prednames, covarnames, fit, vcvres, oneline, version, dependencies, description, vardefs, varunits, filepath) {
    require(yaml)
    list2env(ta.arguments(df, NULL, adm1name, prednames, covarnames, factorouts, demeanfile), environment())

    fp <- file(filepath, 'w')

    variables <- list(outcome=paste0(vardefs[[outname]], " [", varunits[[outname]], "]"))
    for (name in unique(c(prednames, covarnames)))
        variables[[name]] <- paste0(vardefs[[name]], " [", varunits[[name]], "]")

    cat("---\n", fp)
    cat(as.yaml(list(oneline=online, version=version, dependencies=dependencies, description=description,
                     "csvv-version"="girdin-2017-01-10", variables=variables)), fp)
    cat("...\n", fp)

    cat("observations\n")
    cat(nrow(df))
    cat("\nprednames\n")
    cat(paste(prednames, sep=','))
    cat("\ncovarnames\n")
    cat(paste(covarnames, sep=','))

    cat("\ngamma\n")

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
            orderedgammas <- c(orderedgammas, fit$gammas[gammas.so.far + 1])
            vcvorder <- c(vcvorder, length(fit$betas) + gammas.so.far)
        }
    }

    cat(paste(gammas, ','))

    cat("\ngammavcv\n")
    for (ii in length(vcvorder))
        cat(paste(vcvfit$fit[vcvorder[ii], vcvorder], sep=','))

    ## Predict, so that can calculate residual variance
    yyhats <- logspec.predict(xxs, zzs, kls, adm1, fit$betas, fit$gammas, fes)
    cat((yy - yyhat)^2 / (nrow(df) - 2))
}
