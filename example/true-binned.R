##setwd("~/research/gcp/hierarchical-estimation/example")

do.tests <- F

## K = 5 bins currently hard-coded
K.bins <- 5

M.adm1 <- 100 # state regions
NperM.adm2 <- 10 # county regions per state
L.covs <- 2 # covariates (GDP P.C. and Mean T)

Y.years <- 20 # years of observations
D.days <- 365 # days per year

error <- 1

teff0 <- .1 # Size of temperature effect
tstar <- 21 # Lowest-effect temperature

tcenters <- c(12, 16.5, 21, 25.5, 30) # By default, bin 3 has 0 effect (but not necessary)
bineffects <- (tcenters - tstar)^2

## Generate observations
df <- data.frame(rate=c(), bin1=c(), bin2=c(), bin3=c(), bin4=c(), bin5=c(), const=c(), gdppc=c(), meant=c(), adm1=c(), adm2=c())

for (mm in 1:M.adm1) {
    gdppc <- rexp(1, rate=.1) # GDP P.C. in units of $1000 USD
    meant <- rnorm(1, 20, 4)

    for (nn in 1:NperM.adm2) {
        print(c(mm, nn))

        ## Socioeconomic parameters
        rate0 <- rexp(1) # region fixed effect

        ## Climate parameters
        adm2.meant <- meant + rnorm(1, 0, 4)
        sdevt <- rexp(1, .4) + rexp(1, .4)

        df.adm2 <- data.frame()
        for (yy in 1:Y.years) {
            temps <- rnorm(D.days, adm2.meant, sdevt)
            bin1 <- sum(temps < 15)
            bin2 <- sum(temps < 20) - bin1
            bin3 <- sum(temps < 23) - bin2 - bin1
            bin4 <- sum(temps < 27 & temps >= 23)
            bin5 <- sum(temps >= 27)

            rate <- rate0 * exp(-meant / 100 - log(gdppc) / 10) + sum(teff0 * exp(-meant / 200 - log(gdppc) / 20) * bineffects * c(bin1, bin2, bin3, bin4, bin5)) + error * rnorm(1)
            dropeff <- teff0 * bineffects[3] * exp(-meant / 200 - log(gdppc) / 20) * bin3
            df.adm2 <- rbind(df.adm2, data.frame(rate=max(rate, 0), bin1, bin2, bin3, bin4, bin5, dropeff)) # Other columns added below
        }

        df.adm2$adm1 <- mm
        df.adm2$adm2 <- (mm - 1) * NperM.adm2 + nn
        df.adm2$const <- rate0 * exp(-meant / 100 - log(gdppc) / 10)
        df.adm2$gdppc <- gdppc
        df.adm2$meant <- meant

        df <- rbind(df, df.adm2)
    }
}

write.csv(df, "true-binned.csv", row.names=F)

if (do.tests) {
    library(lfe)

    summary(felm(rate ~ 0 + bin1 + bin2 + bin4 + bin5 | adm2 | 0 | adm2, data=df))

    df$loggdppc <- log(df$gdppc)
    summary(felm(rate ~ 0 + bin1 + bin1:loggdppc + bin1:meant + bin2 + bin2:loggdppc + bin2:meant + bin4 + bin4:loggdppc + bin4:meant + bin5 + bin5:loggdppc + bin5:meant | adm2 | 0 | adm2, data=df))
}
