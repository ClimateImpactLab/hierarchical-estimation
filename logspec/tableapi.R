source("logspec.R")
source("search.R")

options(lfe.pint=5)

ta.get.column <- function(df, name) {
    with(df, eval(parse(text=name)))
}

ta.interpret.factorout <- function(df, factorout) {
    if (grepl(':', factorout)) {
        parts <- trimws(strsplit(factorout, ":")[[1]])
        value1 <- ta.get.column(df, parts[1])
        value2 <- ta.get.column(df, parts[2])

        if (is.factor(value1) && is.factor(value2))
            thisfactor <- as.factor(paste(as.character(value1), as.character(value2), sep=':'))
        else if (is.factor(value1)) {
            attr(value1, 'x') <- value2
            thisfactor <- value1
        } else if (is.factor(value2)) {
            attr(value2, 'x') <- value1
            thisfactor <- value2
        } else
            stop(paste("In factorouts, either", parts[1], "or", parts[2], "must be a factor"))
    } else {
        thisfactor <- as.factor(ta.get.column(df, factorout))
    }

    thisfactor
}

ta.arguments.stage1 <- function(df, outname, prednames, factorouts, demeanfile) {
    if (!is.null(outname))
        yy <- df[, outname]

    ## Create the factor list
    if (!is.null(demeanfile)) {
        if (!file.exists(demeanfile))
            stop(paste("Cannot file demean file", demeanfile))
        factors <- demeanfile
    } else if (length(factorouts) == 0) {
        factors <- c()
    } else {
        if (length(factorouts) == 1)
            factors <- ta.interpret.factorout(df, factorouts)
        else {
            factors <- list()
            for (ii in 1:length(factorouts))
                factors[[ii]] <- ta.interpret.factorout(df, factorouts[ii])
        }
    }

    unipreds <- unique(prednames)

    ## Fill in missing predictors
    for (pred in unipreds)
        if (!(pred %in% names(df)))
            df[, pred] <- NA

    xxs <- df[, unipreds]

    if (is.null(outname))
        list(factors=factors, xxs=xxs)
    else
        list(yy=yy, factors=factors, xxs=xxs)
}

ta.arguments <- function(df, outname, adm1name, prednames, covarnames, factorouts, demeanfile) {
    list2env(ta.arguments.stage1(df, outname, prednames, factorouts, demeanfile), environment())

    if (!is.null(adm1name))
        adm1 <- as.numeric(factor(df[, adm1name]))
    else
        adm1 <- NULL

    unipreds <- unique(prednames)
    unicovars <- unique(covarnames)
    if (!('1' %in% unicovars))
        stop("One of the covariates must be '1' for each predictor.")
    unicovars <- unicovars[unicovars != '1'] # Drop this covariate

    ## Generate an order across representative rows
    if (!is.null(adm1name)) {
        zzs.representatives <- !duplicated(adm1)
        zzs.adm1 <- adm1[zzs.representatives]
        zzs.adm1order <- c()
        for (mm in 1:max(adm1))
            zzs.adm1order <- c(zzs.adm1order, which(zzs.adm1 == mm))
        zzs <- df[zzs.representatives, unicovars, drop=F][zzs.adm1order, , drop=F] # pull out, in order of adm1 numbers
    } else
        zzs <- df[, unicovars, drop=F]

    kls <- matrix(F, ncol(xxs), ncol(zzs))
    for (kk in 1:ncol(xxs))
        kls[kk, unicovars %in% covarnames[prednames == unipreds[kk]]] <- T

    if (is.null(outname))
        list(adm1=adm1, factors=factors, xxs=xxs, zzs=zzs, kls=kls)
    else {
        ## Define the standard prior
        zzs.taus <- log(sd(yy)) / apply(zzs, 2, sd)
        taus <- c()
        for (kk in 1:nrow(kls))
            taus <- c(taus, zzs.taus[kls[kk, ]])
        prior <- gaussian.prior(taus)
        gammapriorderiv <- gaussian.gammapriorderiv(taus)

        list(yy=yy, adm1=adm1, factors=factors, xxs=xxs, zzs=zzs, kls=kls, prior=prior, gammapriorderiv=gammapriorderiv)
    }
}

## Wrapper on estimate.logspec
ta.estimate.logspec <- function(df, outname, adm1name, prednames, covarnames, factorouts, weights=1, priorset='default', initset='match', demeanfile=NULL) {
    list2env(ta.arguments(df, outname, adm1name, prednames, covarnames, factorouts, demeanfile), environment())

    ## Variables used by ta.match.marginals (can stay NULL)
    meanzz <- NULL
    meanbetas <- NULL
    mod <- NULL

    ## Pass to search
    known.betas.info <- NULL

    if (priorset == 'none') {
        prior <- noninformative.prior
        gammapriorderiv <- noninformative.gammapriorderiv
    } else if (priorset == 'default' || priorset == 'betaprior') {
        ## Determine mean betas
        mod <- ta.ols(df, outname, adm1name, prednames, covarnames, factorouts, weights=weights)
        meanzz <- get.meanzz(df, covarnames, weights)
        meanbetas <- ta.ols.predict.betas(as.data.frame(t(meanzz)), prednames, covarnames, mod)[1,]

        if (priorset == 'betaprior') {
            require(mvtnorm)

            betaindices <- sapply(unique(prednames), function(pred) which(prednames == pred & covarnames == '1'))
            meanbetavcv <- mod$robustvcv[betaindices, betaindices]
            invmeanbetavcv <- solve(meanbetavcv)

            gaussprior <- prior # prior set by ta.arguments
            gaussgammapriorderiv <- gammapriorderiv

            prior <- function(betas, gammas) {
                obsbetas <- ta.predict.betas(meanzz, prednames, covarnames, betas, gammas)
                dmvnorm(obsbetas, meanbetas, meanbetavcv, log=T) + gaussprior(betas, gammas)
            }

            gammapriorderiv <- function(gammas) {
                ## Need d obj / d gamma = (d obj / d beta) (d beta / d gamma)
                d.dbetas <- invmeanbetavcv %*% (obsbetas - meanbetas) # Kx1
                dbetas.dgammas <- kls * t(matrix(meanzz, length(meanzz), nrow(kls))) # KxL
                t(d.dbetas) %*% dbetas.dgammas + gausspriorderiv(betas, gammas)
            }
        } else {
            ## By default, force betas to match at zz0
            known.betas.info <- list(zz0=meanzz, betas0=meanbetas)
        }
    }

    if (initset == 'match') {
        ## Create initial gammas based on OLS
        print("Finding initial gamma values...")
        initgammas <- ta.match.marginals(df, outname, adm1name, prednames, covarnames, factorouts, zz0=meanzz, betas0=meanbetas, weights, mod.ols=mod, demeanfile=demeanfile)
        print(initgammas)
    } else
        initgammas <- NULL

    search.logspec(yy, xxs, zzs, kls, adm1, factors, weights=weights,
                   initgammas=initgammas, prior=prior,
                   gammapriorderiv=gammapriorderiv, known.betas.info=known.betas.info)
}

## Wrapper on estimate.vcv
ta.estimate.vcv <- function(betas, gammas, sigmas, df, outname, adm1name, prednames, covarnames, factorouts, demeanfile=NULL, ...) {
    list2env(ta.arguments(df, outname, adm1name, prednames, covarnames, factorouts, demeanfile), environment())
    estimate.vcv(betas, gammas, sigmas, yy, xxs, zzs, kls, adm1, factors, prior=prior, ...)
}

ta.predict <- function(df, adm1name, prednames, covarnames, factorouts, betas, gammas, fes=NULL, demeanfile=NULL) {
    list2env(ta.arguments(df, NULL, adm1name, prednames, covarnames, factorouts, demeanfile), environment())

    logspec.predict(xxs, zzs, kls, adm1, betas, gammas, fes)
}

ta.predict.betas <- function(df, prednames, covarnames, betas, gammas) {
    list2env(ta.arguments(df, NULL, NULL, prednames, covarnames, c(), NULL, NULL), environment())

    logspec.predict.betas(zzs, kls, betas, gammas)
}

ta.rsqr <- function(fit, df, outname, adm1name, prednames, covarnames, factorouts, weights=1, demeanfile=NULL) {
    list2env(ta.arguments(df, outname, adm1name, prednames, covarnames, factorouts, demeanfile), environment())

    rsqr(yy, xxs, zzs, kls, adm1, factors, fit$betas, fit$gammas, weights)
}

ta.rsqr.projected <- function(fit, df, outname, adm1name, prednames, covarnames, factorouts, weights=1, demeanfile=NULL) {
    list2env(ta.arguments(df, outname, adm1name, prednames, covarnames, factorouts, demeanfile), environment())

    rsqr.projected(yy, xxs, zzs, kls, adm1, factors, fit$betas, fit$gammas, weights)
}

ta.ols <- function(df, outname, adm1name, prednames, covarnames, factorouts, clustserr=F, weights=1) {
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

    if (clustserr)
        formula <- paste(formula, "|", paste(factorouts, collapse=" + "), "| 0 |", adm1name)
    else
        formula <- paste(formula, "|", paste(factorouts, collapse=" + "))

    print(formula)

    if (length(weights) == 1)
        felm(as.formula(formula), data=df)
    else
        felm(as.formula(formula), weights=weights, data=df)
}

## From https://rdrr.io/github/skranz/regtools/src/R/felm.r
ta.ols.predict <- function(object, newdata, use.fe = TRUE,...) {
    co = coef(object)

    rownames(newdata) = seq_along(newdata[,1])
    ## model matrix without intercept
    mm = model.matrix(object = object,data = newdata)

    if (nrow(mm) < nrow(newdata)) {
        warning("Observations dropped from newdata due to NA.")
    }

    ## remove intercept
    if (ncol(mm) == length(co) + 1) {
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

ta.ols.parameters <- function(prednames, covarnames, mod) {
    betas <- c()
    for (pred in unique(prednames))
        betas <- c(betas, mod$coeff[row.names(mod$coeff) == pred, 1])

    gammaorder <- ta.gammaorder(prednames, covarnames)
    factors <- attributes(mod$terms)$factors > 0

    gammas <- rep(NA, max(gammaorder, na.rm=T))
    for (ii in 1:length(gammaorder)) {
        if (is.na(gammaorder[ii]))
            next

        gammas[gammaorder[ii]] <- mod$coeff[factors[row.names(factors) == prednames[ii],] & factors[row.names(factors) == covarnames[ii],]]
    }

    list(betas=betas, gammas=gammas)
}

ta.ols.predict.betas <- function(df, prednames, covarnames, mod) {
    list2env(ta.arguments(df, NULL, NULL, prednames, covarnames, c(), NULL), environment())

    params <- ta.ols.parameters(prednames, covarnames, mod)
    betas <- params$betas
    gammas <- params$gammas

    result <- matrix(betas, nrow(zzs), length(betas)) # Only modify this for predictors with covariates
    gammas.so.far <- 0 # Keep track of how many coefficients used

    for (kk in 1:nrow(kls)) {
        gammas.here <- sum(kls[kk, ])
        if (gammas.here == 0)
            next # Nothing to do: already dmxxs[, kk]

        mygammas <- gammas[(gammas.so.far+1):(gammas.so.far+gammas.here)]
        gammas.so.far <- gammas.so.far + gammas.here
        result[, kk] <- betas[kk] + as.matrix(zzs[, kls[kk, ]]) %*% mygammas
    }

    result
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
ta.match.marginals <- function(df, outname, adm1name, prednames, covarnames, factorouts, zz0=NULL, betas0=NULL, weights=1, prior=noninformative.prior, mod.ols=NULL, demeanfile=NULL) {
    list2env(ta.arguments(df, outname, adm1name, prednames, covarnames, factorouts, demeanfile), environment())
    list2env(check.arguments(yy, xxs, zzs, kls, adm1, factors), environment())
    list2env(demean.yxs(yy, xxs, factors, weights), environment())

    if (is.null(mod.ols))
        mod.ols <- ta.ols(df, outname, adm1name, prednames, covarnames, factorouts, weights)

    ## Pull out marginals in the right order
    params <- ta.ols.parameters(prednames, covarnames, mod.ols)
    marginals.ols <- params$gammas

    if (is.null(zz0))
        zz0 <- get.meanzz(df, covarnames, weights)
    zzvars <- apply(zzs, 2, sd)

    marginals.zzvars <- c()
    for (kk in 1:nrow(kls))
        marginals.zzvars <- c(marginals.zzvars, zzvars[kls[kk,]])

    objective <- function(gammas) {
        if (is.null(betas0))
            betas0 <- stacked.betas(K, L, gammas, dmyy, dmxxs, zzs, kls, adm1, weights)

        marginals.mle <- marginal.interactions(zz0, kls, betas0, gammas)

        sum((marginals.ols - marginals.mle)^2 / marginals.zzvars)
    }

    ## Decide on initial gammas
    if (is.null(betas0))
        get.betas <- stacked.betas
    else
        get.betas <- make.known.betas(zz0, betas0)

    result <- estimate.logspec.demeaned(dmyy, dmxxs, zzs, kls, adm1, weights, maxiter=1, initgammas=NULL, prior=noninformative.prior, get.betas=get.betas)

    opt <- optim(result$gammas, objective)
    if (opt$convergence == 0)
        return(NULL)

    gammas <- opt$par

    # Compare to all 0's
    betas <- get.betas(K, L, gammas * 0, dmyy, dmxxs, zzs, kls, adm1, weights)
    likeli0 <- calc.likeli.demeaned(dmxxs, dmyy, zzs, kls, adm1, betas, gammas * 0,
                                    result$sigmas, weights, prior)

    ## Decrease if too extreme
    for (ii in 1:10) { # To 0.1% of gammas
        betas <- get.betas(K, L, gammas, dmyy, dmxxs, zzs, kls, adm1, weights)
        likeli <- calc.likeli.demeaned(dmxxs, dmyy, zzs, kls, adm1, betas, gammas, result$sigmas, weights, prior)
        if (likeli > likeli0) {
            print(ii)
            return(gammas)
        }
        gammas <- gammas / 2
    }

    gammas * 0
}

get.meanzz <- function(df, covarnames, weights=1) {
    meanzz <- c()
    zzcovarnames <- c()

    for (covar in unique(covarnames)) {
        if (covar != '1') {
            zzcovarnames <- c(zzcovarnames, covar)
            if (length(weights) == 1)
                meanzz <- c(meanzz, mean(df[, covar]))
            else
                meanzz <- c(meanzz, weighted.mean(df[, covar], weights))
        }
    }

    names(meanzz) <- zzcovarnames
    meanzz
}

ta.save.demeaned <- function(filename, df, outname, prednames, factorouts, weights=1) {
    list2env(ta.arguments.stage1(df, outname, prednames, factorouts, demeanfile=NULL), environment())

    save.demeaned(yy, xxs, factors, weights, filename)
}
