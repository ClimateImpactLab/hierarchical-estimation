setwd("~/research/gcp/hierarchical-estimation/rtree")
source("../../socioeconomics/helpers/header.R")

alltbl <- NULL
for (filename in c("bin_13C_18C.csv", "bin_23C_28C.csv", "bin_28C_33C.csv", "bin_33C_InfC.csv", "bin_3C_8C.csv", "bin_8C_13C.csv", "bin_n12C_n7C.csv", "bin_n17C_n12C.csv", "bin_n2C_3C.csv", "bin_n7C_n2C.csv", "bin_nInfC_n17C.csv")) {
    tbl <- read.csv(header.deparse(paste0("../../data/adaptation/predictors-space-all/", filename)))
    names(tbl)[3:4] <- paste0(filename, c('_coef', '_serr'))
    ##tbl[is.nan(tbl[, 3]), 3] <- 0
    if (is.null(alltbl)) {
        alltbl <- tbl
    } else {
        alltbl[, paste0(filename, c('_coef', '_serr'))] <- tbl[, paste0(filename, c('_coef', '_serr'))]
    }
}

alltbl$loggdppc <- log(alltbl$gdppc)
alltbl$logpopop <- log(alltbl$popop)
1
library(party)

## weights <- (1*!is.na(alltbl$bin_13C_18C.csv_coef))/alltbl$bin_13C_18C.csv_serr^2 + (1*!is.na(alltbl$bin_23C_28C.csv_coef))/alltbl$bin_23C_28C.csv_serr^2 + (1*!is.na(alltbl$bin_28C_33C.csv_coef))/alltbl$bin_28C_33C.csv_serr^2 + (1*!is.na(alltbl$bin_33C_InfC.csv_coef))/alltbl$bin_33C_InfC.csv_serr^2 + (1*!is.na(alltbl$bin_3C_8C.csv_coef))/alltbl$bin_3C_8C.csv_serr^2 + (1*!is.na(alltbl$bin_8C_13C.csv_coef))/alltbl$bin_8C_13C.csv_serr^2 + (1*!is.na(alltbl$bin_n12C_n7C.csv_coef))/alltbl$bin_n12C_n7C.csv_serr^2 + (1*!is.na(alltbl$bin_n17C_n12C.csv_coef))/alltbl$bin_n17C_n12C.csv_serr^2 + (1*!is.na(alltbl$bin_n2C_3C.csv_coef))/alltbl$bin_n2C_3C.csv_serr^2 + (1*!is.na(alltbl$bin_n7C_n2C.csv_coef))/alltbl$bin_n7C_n2C.csv_serr^2 + (1*!is.na(alltbl$bin_nInfC_n17C.csv_coef))/alltbl$bin_nInfC_n17C.csv_serr^2

## rt <- ctree(bin_13C_18C.csv_coef + bin_23C_28C.csv_coef + bin_28C_33C.csv_coef + bin_33C_InfC.csv_coef + bin_3C_8C.csv_coef + bin_8C_13C.csv_coef + bin_n12C_n7C.csv_coef + bin_n17C_n12C.csv_coef + bin_n2C_3C.csv_coef + bin_n7C_n2C.csv_coef + bin_nInfC_n17C.csv_coef ~ group + loggdppc + logpopop + meandays_nInfC_n17C + meandays_n17C_n12C + meandays_n12C_n7C + meandays_n7C_n2C + meandays_n2C_3C + meandays_3C_8C + meandays_8C_13C + meandays_13C_18C + meandays_23C_28C + meandays_28C_33C + meandays_33C_InfC, data=alltbl[valid,], weights=weights[valid], maxdepth=3)

## plot(rt)

## valid <- !is.na(alltbl$bin_28C_33C.csv_coef) & !is.na(alltbl$bin_33C_InfC.csv_coef)

## ct <- ctree(bin_28C_33C.csv_coef + bin_33C_InfC.csv_coef ~ group + loggdppc + logpopop + meandays_nInfC_n17C + meandays_n17C_n12C + meandays_n12C_n7C + meandays_n7C_n2C + meandays_n2C_3C + meandays_3C_8C + meandays_8C_13C + meandays_13C_18C + meandays_23C_28C + meandays_28C_33C + meandays_33C_InfC, data=alltbl[valid,], controls= ctree_control(mincriterion = 0.5))

## plot(ct, type="simple")

library(rpart)
library(rpart.plot)

fit <- rpart(bin_33C_InfC.csv_coef ~ group + loggdppc + logpopop + meantemp + meandays_33C_InfC, data=alltbl, weights=1/alltbl$bin_33C_InfC.csv_serr^2)

prp(fit, type=4, extra=1, yesno=T, fallen.leaves=T, faclen=0, varlen=0, branch.lty=3)

fit <- rpart(bin_28C_33C.csv_coef ~ group + loggdppc + logpopop + meantemp + meandays_28C_33C, data=alltbl, weights=1/alltbl$bin_28C_33C.csv_serr^2)

prp(fit, type=4, extra=1, yesno=T, fallen.leaves=T, faclen=0, varlen=0, branch.lty=3)

fit <- rpart(bin_n17C_n12C.csv_coef ~ group + loggdppc + logpopop + meantemp + meandays_n17C_n12C, data=alltbl, weights=1/alltbl$bin_n17C_n12C.csv_serr^2)

prp(fit, type=4, extra=1, yesno=T, fallen.leaves=T, faclen=0, varlen=0, branch.lty=3)
