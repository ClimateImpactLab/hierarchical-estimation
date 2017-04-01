##setwd("~/research/gcp/hierarchical-estimation/example")

M.adm1 <- 100 # state regions
NperM.adm2 <- 10 # county regions per state

Y.years <- 20 # years of observations
D.days <- 365 # days per year

teff0 <- .1 # Size of temperature effect
tstar <- 22 # Lowest-effect temperature

## Generate observations
df <- data.frame(rate=c(), sumlev=c(), sumsqr=c(), gdppc=c(), meant=c(), adm1=c(), adm2=c())

for (mm in 1:M.adm1) {
    for (nn in 1:NperM.adm2) {
        print(c(mm, nn))

        ## Socioeconomic parameters
        rate0 <- rexp(1) # region fixed effect

        ## Climate parameters
        meant <- rnorm(1, 20, 5)
        sdevt <- rexp(1, .4) + rexp(1, .4)

        gdppc <- rexp(1, rate=.1 + abs(meant - tstar) * sdevt / 10) # GDP P.C. in units of $1000 USD

        df.adm2 <- data.frame()
        for (yy in 1:Y.years) {
            temps <- rnorm(D.days, meant, sdevt)
            sumlev <- sum(temps)
            sumsqr <- sum(temps^2)

            rate <- rate0 * exp(-log(gdppc) / 10 - meant / 100) + sum(teff0 * exp(-log(gdppc) / 20 - meant / 200) * (temps - tstar)^2) + rnorm(1)
            df.adm2 <- rbind(df.adm2, data.frame(rate=max(rate, 0), sumlev, sumsqr)) # Other columns added below
        }

        df.adm2$adm1 <- mm
        df.adm2$adm2 <- (mm - 1) * NperM.adm2 + nn
        df.adm2$gdppc <- gdppc
        df.adm2$meant <- meant

        df <- rbind(df, df.adm2)
    }
}

library(lfe)

summary(felm(rate ~ 0 + sumlev + sumsqr | adm2 | 0 | adm2, data=df))

df$loggdppc <- log(df$gdppc)
mod <- felm(rate ~ 0 + sumlev + sumlev:loggdppc + sumlev:meant + sumsqr + sumsqr:loggdppc + sumsqr:meant | adm2 | 0 | adm2, data=df)
summary(mod)

write.csv(df, "quadratic.csv", row.names=F)

makecurve <- function(mod, loggdppc, meant, tseq) {
    beta.lev <- mod$coeff[1] + mod$coeff[3] * loggdppc + mod$coeff[4] * meant
    beta.sqr <- mod$coeff[2] + mod$coeff[5] * loggdppc + mod$coeff[6] * meant
    beta.lev * (tseq - tstar) + beta.sqr * (tseq^2 - tstar^2)
}

tseq <- seq(10, 35, length.out=100)
plotdf <- data.frame(temp=rep(tseq, 5),
                     resp=c(makecurve(mod, quantile(df$loggdppc, .1), quantile(df$meant, .1), tseq),
                            makecurve(mod, quantile(df$loggdppc, .1), quantile(df$meant, .9), tseq),
                            makecurve(mod, quantile(df$loggdppc, .5), quantile(df$meant, .5), tseq),
                            makecurve(mod, quantile(df$loggdppc, .9), quantile(df$meant, .1), tseq),
                            makecurve(mod, quantile(df$loggdppc, .9), quantile(df$meant, .9), tseq)),
                     group=rep(c("low-low", "low-high", "median", "high-low", "high-high"), each=100))

library(ggplot2)

ggplot(plotdf, aes(temp, resp, colour=group)) +
    geom_line()
