## K = 5 bins currently hard-coded
K.bins <- 5

M.adm1 <- 100 # state regions
NperM.adm2 <- 10 # county regions per state
L.covs <- 1 # covariates (just GDP P.C.)

Y.years <- 20 # years of observations
D.days <- 365 # days per year

teff0 <- .1 # Size of temperature effect
tstar <- 22 # Lowest-effect temperature

## Generate observations
df <- data.frame(rate=c(), bin1=c(), bin2=c(), bin3=c(), bin4=c(), bin5=c(), gdppc=c(), meant=c(), adm1=c(), adm2=c())

for (mm in 1:M.adm1) {
    for (nn in 1:NperM.adm2) {
        print(c(mm, nn))

        ## Socioeconomic parameters
        rate0 <- rexp(1) # region fixed effect
        gdppc <- rexp(1, rate=.1) # GDP P.C. in units of $1000 USD

        ## Climate parameters
        meant <- rnorm(1, 20, 5)
        sdevt <- rexp(1, .4) + rexp(1, .4)

        df.adm2 <- data.frame()
        for (yy in 1:Y.years) {
            temps <- rnorm(D.days, meant, sdevt)
            bin1 <- sum(temps < 15)
            bin2 <- sum(temps < 20) - bin1
            bin3 <- sum(temps < 23) - bin2 - bin1
            bin4 <- sum(temps < 27 & temps >= 23)
            bin5 <- sum(temps >= 27)

            rate <- rate0 * exp(-gdppc / 10 - meant / 100) + sum(teff0 * exp(-gdppc / 20 - meant / 200) * (temps - tstar)^2) + rnorm(1)
            df.adm2 <- rbind(df.adm2, data.frame(rate=max(rate, 0), bin1, bin2, bin3, bin4, bin5)) # Other columns added below
        }

        df.adm2$adm1 <- mm
        df.adm2$adm2 <- (mm - 1) * NperM.adm2 + nn
        df.adm2$gdppc <- gdppc
        df.adm2$meant <- meant

        df <- rbind(df, df.adm2)
    }
}

library(lfe)

summary(felm(rate ~ 0 + bin1 + bin2 + bin4 + bin5 | adm2 | 0 | adm2, data=df))

summary(felm(rate ~ 0 + bin1 + bin1:gdppc + bin1:meant + bin2 + bin2:gdppc + bin2:meant + bin4 + bin4:gdppc + bin4:meant + bin5 + bin5:gdppc + bin5:meant | adm2 | 0 | adm2, data=df))

write.csv(df, "binned.csv", row.names=F)
