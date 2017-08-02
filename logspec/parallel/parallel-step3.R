##setwd("~/research/gcp/hierarchical-estimation/interpolate")

## Load the library
source("logspec.R")

## Load the data
df <- read.csv("../example/true-binned.csv")
df$log_gdppc <- log(df$gdppc)

methast.result <- parallel.combine.methast("example", 10,
                                           df$rate, df[, c('bin1', 'bin2', 'bin4', 'bin5')],
                                           df[!duplicated(df$adm1), c('meant', 'log_gdppc', 'meant', 'log_gdppc',
                                                                      'meant', 'log_gdppc', 'meant', 'log_gdppc')],
                                           df$adm1, df$adm2, prior) # need to define prior

save(methast.result, file=paste0("MH-combined.RData"))
