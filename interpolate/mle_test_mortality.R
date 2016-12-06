setwd("/Users/trinettachong/Dropbox/MLE/")
#clear environment
rm(list = ls())

#Load packages
library(readstata13)
library(RcppArmadillo)
library(nnls)
library(lfe)

#Load data (mortality with 6 countries)
source("logspec.R")
df <- read.dta13("mle_mortality_5dec.dta")
df$log_gdppc <- log(df$gdppc)
df$log_popop <- log(df$popop)

#view all variables in dataframe
names(df)

#create a vector for each bin
df$bin1 <- df$bin_nInfC_n17C_BEST 
df$bin2 <- df$bin_n17C_n12C_BEST 
df$bin3 <- df$bin_n12C_n7C_BEST 
df$bin4 <- df$bin_n7C_n2C_BEST 
df$bin5 <- df$bin_n2C_3C_BEST 
df$bin6 <- df$bin_3C_8C_BEST 
df$bin7 <- df$bin_8C_13C_BEST 
df$bin8 <- df$bin_13C_18C_BEST 
df$bin9 <- df$bin_18C_23C_BEST 
df$bin10 <- df$bin_23C_28C_BEST 
df$bin11 <- df$bin_28C_33C_BEST 
df$bin12 <- df$bin_33C_InfC_BEST

#check that bin is a vector
is.vector(df$bin10)

## Find the coefficients (2 covariates, 11 bins)
#MLE
timetoresult <- proc.time()
result <- estimate.logspec(df$deathrate0, df[, c('bin1', 'bin2', 'bin3', 'bin4', 'bin5', 'bin6', 'bin7', 'bin8', 'bin10', 'bin11', 'bin12' )],
                           df[!duplicated(df$adm1_mle), c('log_popop', 'log_gdppc', 'log_popop', 'log_gdppc',
                                                          'log_popop', 'log_gdppc', 'log_popop', 'log_gdppc',
                                                          'log_popop', 'log_gdppc', 'log_popop', 'log_gdppc',
                                                          'log_popop', 'log_gdppc', 'log_popop', 'log_gdppc',
                                                          'log_popop', 'log_gdppc', 'log_popop', 'log_gdppc',
                                                          'log_popop', 'log_gdppc')],
                           df$adm1_mle, df$adm2_mle)
proc.time() - timetoresult
print(result)

#OLS
result_ols <- felm(deathrate0 ~
                     bin1 + bin1:log_popop + bin1:log_gdppc 
                   + bin2 + bin2:log_popop + bin2:log_gdppc 
                   + bin3 + bin3:log_popop + bin3:log_gdppc
                   + bin4 + bin4:log_popop + bin4:log_gdppc 
                   + bin5 + bin5:log_popop + bin5:log_gdppc 
                   + bin6 + bin6:log_popop + bin6:log_gdppc 
                   + bin7 + bin7:log_popop + bin7:log_gdppc 
                   + bin8 + bin8:log_popop + bin8:log_gdppc 
                   + bin10 + bin10:log_popop + bin10:log_gdppc 
                   + bin11 + bin11:log_popop + bin11:log_gdppc 
                   + bin12 + bin12:log_popop + bin12:log_gdppc| adm2_mle | 0 | adm2_mle, data=df)
summary(result_ols)


## Find the coefficients (2 covariates, 2 bins)
#MLE
timetoresult <- proc.time()
result <- estimate.logspec(df$deathrate0, df[, c('bin10', 'bin11')],
                           df[!duplicated(df$adm1_mle), c('log_popop', 'log_gdppc', 'log_popop', 'log_gdppc')],
                           df$adm1_mle, df$adm2_mle)
proc.time() - timetoresult
print(result)

#OLS
result_ols <- felm(deathrate0 ~ 0 + bin10 + bin10:log_popop + bin10:log_gdppc + bin11 + bin11:log_popop + bin11:log_gdppc | adm2_mle | 0 | adm2_mle, data=df)
summary(result_ols)
