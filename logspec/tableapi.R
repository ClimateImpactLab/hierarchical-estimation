source("logspec.R")
source("search.R")

ta.arguments <- function(df, outname, adm1name, adm2name, prednames, covarnames) {
    if (!is.null(outname))
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
    zzs.adm1 <- adm1[zzs.representatives]
    zzs.adm1order <- c()
    for (mm in 1:max(adm1))
        zzs.adm1order <- c(zzs.adm1order, which(zzs.adm1 == mm))
    zzs <- df[zzs.representatives, unicovars, drop=F][zzs.adm1order, , drop=F] # pull out, in order of adm1 numbers

    kls <- matrix(F, ncol(xxs), ncol(zzs))
    for (kk in 1:ncol(xxs))
        kls[kk, unicovars %in% covarnames[prednames == unipreds[kk]]] <- T

    if (is.null(outname))
        list(adm1=adm1, adm2=adm2, xxs=xxs, zzs=zzs, kls=kls)
    else {
        ## Define the standard prior
        zzs.taus <- log(sd(yy)) / apply(zzs, 2, sd)
        taus <- c()
        for (kk in 1:nrow(kls))
            taus <- c(taus, zzs.taus[kls[kk, ]])
        gammaprior <- gaussian.gammaprior(taus)

        list(yy=yy, adm1=adm1, adm2=adm2, xxs=xxs, zzs=zzs, kls=kls, gammaprior=gammaprior)
    }
}

## Wrapper on estimate.logspec
ta.estimate.logspec <- function(df, outname, adm1name, adm2name, prednames, covarnames, weights=1, priorset='default') {
    list2env(ta.arguments(df, outname, adm1name, adm2name, prednames, covarnames), environment())

    if (priorset == 'none')
        gammaprior <- noninformative.gammaprior

    ## Create initial gammas based on OLS
    print("Finding initial gamma values...")
    initgammas <- ta.match.marginals(df, outname, adm1name, adm2name, prednames, covarnames, weights)

    search.logspec(yy, xxs, zzs, kls, adm1, adm2, weights=weights,
                   initgammas=initgammas, gammaprior=gammaprior)
}

## Wrapper on estimate.vcv
ta.estimate.vcv <- function(betas, gammas, sigmas, df, outname, adm1name, adm2name, prednames, covarnames, ...) {
    list2env(ta.arguments(df, outname, adm1name, adm2name, prednames, covarnames), environment())
    estimate.vcv(betas, gammas, sigmas, yy, xxs, zzs, kls, adm1, adm2, gammaprior=gammaprior, ...)
}

ta.predict <- function(df, adm1name, adm2name, prednames, covarnames, betas, gammas, fes=NULL) {
    list2env(ta.arguments(df, NULL, adm1name, adm2name, prednames, covarnames), environment())

    logspec.predict(xxs, zzs, kls, adm1, adm2, betas, gammas, fes)
}

ta.rsqr <- function(fit, df, outname, adm1name, adm2name, prednames, covarnames, weights=1) {
    list2env(ta.arguments(df, outname, adm1name, adm2name, prednames, covarnames), environment())

    rsqr(yy, xxs, zzs, kls, adm1, adm2, fit$betas, fit$gammas, weights)
}

ta.rsqr.projected <- function(fit, df, outname, adm1name, adm2name, prednames, covarnames, weights=1) {
    list2env(ta.arguments(df, outname, adm1name, adm2name, prednames, covarnames), environment())

    rsqr.projected(yy, xxs, zzs, kls, adm1, adm2, fit$betas, fit$gammas, weights)
}

ta.ols <- function(df, outname, adm1name, adm2name, prednames, covarnames, weights=1) {
    require(lfe)
    if (length(prednames) != length(covarnames))
        stop("prednames and covarnames must be the same length.")

    formula <- paste(outname, "~ 0")
    for (ii in 1:length(prednames)) {
        if (covarnames[ii] == '1')
            formula <- paste(formula, "+", prednames[ii])
        else
            formula <- paste(formula, "+", paste(prednames[ii], covarnames[ii], sep=':'))
    }

    formula <- paste(formula, "|", adm2name, "| 0 |", adm1name)

    felm(as.formula(formula), data=df)
}

## From https://rdrr.io/github/skranz/regtools/src/R/felm.r
ta.ols.predict <- function(object, newdata, use.fe = TRUE,...) {
    co = coef(object)

    rownames(newdata) = seq_along(newdata[,1])
    ## model matrix without intercept
    mm = model.matrix(object = object,data = newdata)

    if (NROW(mm)<NROW(newdata)) {
        warning("Observations dropped from newdata due to NA.")
    }

    ## remove intercept
    if (NCOL(mm)==length(co)+1) {
        mm = mm[,-1,drop=FALSE]
    }
    y.pred = mm %*% co

    fe.vars = names(object$fe)
    if (use.fe & length(fe.vars)>0) {
        rows = as.integer(rownames(mm))
        nd = newdata[rows,]
        all.fe = getfe(object)
        fe.var = fe.vars[1]
        for (fe.var in fe.vars) {
            df = all.fe[all.fe$fe == fe.var,]
            frows = match(nd[[fe.var]],df$idx)
            myfe = df$effect[frows]
            myfe[is.na(myfe)] = 0

            y.pred = y.pred + myfe
        }
    }
    as.vector(y.pred)
}


## Return the set of gamma indices that correspond to each predname-covarname combination
ta.gammaorder <- function(prednames, covarnames) {
    if (length(prednames) != length(covarnames))
        stop("prednames and covarnames must be the same length.")

    unipreds <- unique(prednames)
    unicovars <- unique(covarnames)
    if (!('1' %in% unicovars))
        stop("One of the covariates must be '1' for each predictor.")
    unicovars <- unicovars[unicovars != '1'] # Drop this covariate

    kls <- matrix(F, length(unipreds), length(unicovars))
    for (kk in 1:length(unipreds))
        kls[kk, unicovars %in% covarnames[prednames == unipreds[kk]]] <- T

    kls.index <- matrix(0, nrow(kls), ncol(kls))
    gammas.so.far <- 0
    for (kk in 1:nrow(kls)) {
        gammas.here <- sum(kls[kk, ])
        if (gammas.here == 0)
            next
        kls.index[kk, kls[kk, ]] <- (gammas.so.far+1):(gammas.so.far+gammas.here)

        gammas.so.far <- gammas.so.far + gammas.here
    }

    indices <- rep(NA, length(prednames))
    for (ii in 1:length(prednames))
        if (covarnames[ii] != '1')
            indices[ii] <- kls.index[which(unipreds == prednames[ii]), which(unicovars == covarnames[ii])]

    indices
}

## Find the set of gammas that has the same marginal effects
ta.match.marginals <- function(df, outname, adm1name, adm2name, prednames, covarnames, weights=1) {
    list2env(ta.arguments(df, outname, adm1name, adm2name, prednames, covarnames), environment())
    list2env(check.arguments(yy, xxs, zzs, kls, adm1, adm2), environment())
    list2env(demean.yxs(K, yy, xxs, adm2), environment())

    mod.ols <- ta.ols(df, outname, adm1name, adm2name, prednames, covarnames, weights)
    ## Pull out marginals in the right order
    marginals.ols <- rep(NA, sum(covarnames != '1'))
    gammaorder <- ta.gammaorder(prednames, covarnames)
    for (ii in 1:length(gammaorder))
        if (!is.na(gammaorder[ii]))
            marginals.ols[gammaorder[ii]] <- mod.ols$coeff[ii]

    zzmeans <- colMeans(zzs)
    zzvars <- apply(zzs, 2, sd)

    marginals.zzvars <- c()
    for (kk in 1:nrow(kls))
        marginals.zzvars <- c(marginals.zzvars, zzvars[kls[kk,]])

    objective <- function(gammas) {
        betas <- stacked.betas(K, L, gammas, dmyy, dmxxs, zzs, kls, adm1, weights)

        marginals.mle <- marginal.interactions(zzmeans, kls, betas, gammas)

        sum((marginals.ols - marginals.mle)^2 / marginals.zzvars)
    }

    ## Decide on initial gammas
    result <- estimate.logspec.demeaned(dmyy, dmxxs, zzs, kls, adm1, adm2, weights, maxiter=1, initgammas=NULL, gammaprior=noninformative.gammaprior)

    optim(result$gammas, objective)$par
}

