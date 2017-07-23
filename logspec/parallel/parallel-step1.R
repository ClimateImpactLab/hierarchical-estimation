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
                           df$adm1, df$adm2, gammaprior=gammaprior)

ses <- estimate.se(result$betas, result$gammas, result$sigmas,
                   df$rate, df[, c('bin1', 'bin2', 'bin4', 'bin5')],
                   df[!duplicated(df$adm1), c('meant', 'log_gdppc', 'meant', 'log_gdppc',
                                              'meant', 'log_gdppc', 'meant', 'log_gdppc')],
                   df$adm1, df$adm2, gammaprior=gammaprior) # define gammaprior

result$betaerr <- ses[1:length(result$betas)]
result$gammaerr <- matrix(ses[(length(result$betas)+1):(length(result$betas)+length(result$gammas))], length(result$betas), length(result$gammas) / length(result$betas))

save(result, file="example-parallel.RData")
