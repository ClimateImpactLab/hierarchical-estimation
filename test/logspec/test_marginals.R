library(RUnit)

## Load the library
source("logspec/tableapi.R", chdir=T)

test.marginal.match <- function() {
    df <- read.csv("example/true-binned.csv")
    df$log_gdppc <- log(df$gdppc)

    prednames <- c('bin1', 'bin1', 'bin1', 'bin2', 'bin2', 'bin2',
                   'bin4', 'bin4', 'bin4', 'bin5', 'bin5', 'bin5')
    covarnames <- c('1', 'meant', 'log_gdppc', '1', 'meant', 'log_gdppc',
                    '1', 'meant', 'log_gdppc', '1', 'meant', 'log_gdppc')

    list2env(ta.arguments(df, 'rate', 'adm1', 'adm2', prednames, covarnames), environment())

    mod <- ta.ols(df, 'rate', 'adm1', 'adm2', prednames, covarnames)
    meanzz <- get.meanzz(df, covarnames, zzs)
    meanbetas <- ta.ols.predict.betas(as.data.frame(t(meanzz)), prednames, covarnames, mod)[1,]

    gammas <- ta.match.marginals(df, 'rate', 'adm1', 'adm2', prednames, covarnames,
                                 zz0=meanzz, betas0=meanbetas)

    betas <- make.known.betas(meanzz, meanbetas)(K, L, gammas, dmyy, dmxxs.orig, zzs, kls, mm, 1)

    betas.mle <- ta.predict.betas(as.data.frame(t(meanzz)), prednames, covarnames, betas, gammas)
    betas.ols <- ta.ols.predict.betas(as.data.frame(t(meanzz)), prednames, covarnames, mod)

    checkEqualsNumeric(betas.mle, betas.ols)
}

