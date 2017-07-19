library(MASS)
library(xtable)

## Don't include '1' in mycovarnames
read.csvv <- function(filepath, myprednames, mycovarnames, covars, evalpts, mciters=1000, probs=c(.025, .5, .975)) {
    fp <- file(filepath, "r")
    while (T) {
        line = readLines(fp, n=1)
        if (length(line) == 0)
            break

        if (line == 'prednames')
            prednames <- strsplit(readLines(fp, n=1), ',')[[1]]
        if (line == 'covarnames')
            covarnames <- strsplit(readLines(fp, n=1), ',')[[1]]
        if (line == 'gamma')
            gamma <- as.numeric(strsplit(readLines(fp, n=1), ',')[[1]])
        if (line == 'gammavcv')
            gammavcv <- t(sapply(strsplit(readLines(fp, n=length(gamma)), ','), as.numeric))
    }
    close(fp)

    ## Generate MC draws from beta(covars)
    if (mciters == 1)
        gammas <- matrix(gamma, 1, length(gamma)) # Give median
    else
        gammas <- mvrnorm(mciters, gamma, gammavcv)
    betas <- data.frame(iter=1:mciters)
    for (predname in myprednames) {
        predtotal <- gammas[, prednames == predname & covarnames == '1']
        for (covarname in mycovarnames)
            predtotal <- predtotal + gammas[, prednames == predname & covarnames == covarname]
        betas[, predname] <- predtotal
    }
    betas <- as.matrix(betas[,-1])

    ## Determine the confidence intervals
    allpoints <- evalpts(betas[1,])
    if (mciters > 1) {
        for (ii in 2:mciters)
            allpoints <- rbind(allpoints, evalpts(betas[ii,]))
    }

    if (mciters == 1)
        allpoints
    else
        apply(allpoints, 2, function(col) quantile(col, probs=probs))
}

latex.csvv <- function(filepath, digits=4) {
    fp <- file(filepath, "r")
    while (T) {
        line = readLines(fp, n=1)
        if (length(line) == 0)
            break

        if (line == 'prednames')
            prednames <- strsplit(readLines(fp, n=1), ',')[[1]]
        if (line == 'covarnames')
            covarnames <- strsplit(readLines(fp, n=1), ',')[[1]]
        if (line == 'gamma')
            gamma <- as.numeric(strsplit(readLines(fp, n=1), ',')[[1]])
        if (line == 'gammavcv')
            gammavcv <- t(sapply(strsplit(readLines(fp, n=length(gamma)), ','), as.numeric))
    }
    close(fp)

    df <- data.frame(variable=c(), gamma=c(), stderr=c())
    for (predname in unique(prednames))
        for (covarname in unique(covarnames)) {
            if (covarname == '1')
                variable <- predname
            else
                variable <- paste(predname, covarname, sep=':')

            kk <- which(prednames == predname & covarnames == covarname)
            df <- rbind(df, data.frame(variable, gamma[kk], paste0('(', format(sqrt(gammavcv[kk, kk]), digits=digits), ')')))
        }

    names(df) <- c("Variable", "Gamma", "Std. Err.")
    xtable(df, digits=digits)
}
