##setwd("~/research/gcp/hierarchical-estimation/interpolate")

## Load the library
source("logspec.R")

## Load the data
df <- read.csv("../example/true-binned.csv")
df$log_gdppc <- log(df$gdppc)

## Find the coefficients
result <- estimate.logspec(df$rate, df[, c('bin1', 'bin2', 'bin4', 'bin5')],
                           df[!duplicated(df$adm1), c('meant', 'log_gdppc', 'meant', 'log_gdppc',
                                                      'meant', 'log_gdppc', 'meant', 'log_gdppc')],
                           df$adm1, df$adm2)

## Try to solve without hierarchy
estimate.logspec(df$rate, df[, c('bin1', 'bin2', 'bin4', 'bin5')],
                 df[!duplicated(df$adm2), c('meant', 'log_gdppc', 'meant', 'log_gdppc',
                                            'meant', 'log_gdppc', 'meant', 'log_gdppc')],
                 df$adm2, df$adm2)

## Solve using R optimization method
estimate.logspec.optim(df$rate, df[, c('bin1', 'bin2', 'bin4', 'bin5')],
                       df[!duplicated(df$adm2), c('meant', 'log_gdppc', 'meant', 'log_gdppc',
                                                  'meant', 'log_gdppc', 'meant', 'log_gdppc')],
                       df$adm2, df$adm2)
