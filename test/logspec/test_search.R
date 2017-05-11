library(RUnit)

## Load the library
source("interpolate/search.R", chdir=T)

teff0 <- .1 # Size of temperature effect
tstar <- 21 # Lowest-effect temperature

tcenters <- c(12, 16.5, 21, 25.5, 30) # By default, bin 3 has 0 effect (but not necessary)
bineffects <- teff0 * (tcenters - tstar)^2


test.estimate.logspec <- function() {
    ## Load the data
    df <- read.csv("example/true-binned.csv")
    df$log_gdppc <- log(df$gdppc)

    result <- search.logspec(df$rate, df[, c('bin1', 'bin2', 'bin4', 'bin5')],
                             df[!duplicated(df$adm1), c('meant', 'log_gdppc')],
                             matrix(T, 4, 2), df$adm1, df$adm2)

    checkEqualsNumeric(result$betas, bineffects[-3], tolerance=.01)

    gamma1 <- -1/200
    gamma2 <- -1/20
    checkEqualsNumeric(result$gammas, rep(c(gamma1, gamma2), 4), tolerance=.01)
}
