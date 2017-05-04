##setwd("~/research/gcp/hierarchical-estimation/logspec")

subdf <- read.csv("mortality_sub5000.csv")

## Prepare the data
subdf$log_gdppc <- log(subdf$gdppc)
subdf$log_popop <- log(subdf$popop)
subdf <- subdf[rowSums(is.na(subdf[, c('dr_age', 'bin_nInf_n17C_pop', 'bin_n17C_n12C_pop', 'bin_n12C_n7C_pop', 'bin_n7C_n2C_pop', 'bin_n2C_3C_pop', 'bin_3C_8C_pop', 'bin_8C_13C_pop', 'bin_13C_18C_pop', 'bin_23C_28C_pop', 'bin_28C_33C_pop', 'bin_33C_Inf_pop', 'log_gdppc', 'log_popop')])) == 0,]

subdf$adm1 <- factor(subdf$adm1)
subdf$adm2 <- factor(subdf$adm2)

stan.data <- list(I = nrow(subdf), N = length(unique(subdf$adm2)),
                  M = length(unique(subdf$adm1)), K = 11, L = 3,
                  nn = as.numeric(subdf$adm2), mm = as.numeric(subdf$adm1),
                  xx = t(subdf[, c('bin_nInf_n17C_pop', 'bin_n17C_n12C_pop', 'bin_n12C_n7C_pop', 'bin_n7C_n2C_pop', 'bin_n2C_3C_pop', 'bin_3C_8C_pop', 'bin_8C_13C_pop', 'bin_13C_18C_pop', 'bin_23C_28C_pop', 'bin_28C_33C_pop', 'bin_33C_Inf_pop')]),
                  zz = aperm(array(unlist(list(subdf[, c('meandays_nInf_n17C', 'log_gdppc', 'log_popop')],
                                         subdf[, c('meandays_n17C_n12C', 'log_gdppc', 'log_popop')],
                                         subdf[, c('meandays_n12C_n7C', 'log_gdppc', 'log_popop')],
                                         subdf[, c('meandays_n7C_n2C', 'log_gdppc', 'log_popop')],
                                         subdf[, c('meandays_n2C_3C', 'log_gdppc', 'log_popop')],
                                         subdf[, c('meandays_3C_8C', 'log_gdppc', 'log_popop')],
                                         subdf[, c('meandays_8C_13C', 'log_gdppc', 'log_popop')],
                                         subdf[, c('meandays_13C_18C', 'log_gdppc', 'log_popop')],
                                         subdf[, c('meandays_23C_28C', 'log_gdppc', 'log_popop')],
                                         subdf[, c('meandays_28C_33C', 'log_gdppc', 'log_popop')],
                                         subdf[, c('meandays_33C_InfC', 'log_gdppc', 'log_popop')])),
                                   dim = c(nrow(subdf), 3, 11)), c(3, 1, 2)),
                  yy = subdf$dr_age)
