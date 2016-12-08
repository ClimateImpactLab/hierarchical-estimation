##setwd("~/research/gcp/hierarchical-estimation/interpolate")

source("../interpolate/logspec.R")

result <- methast(500, .1, .1, function(v) dnorm(v, log=T))

result <- repeated.methast(5, 600, 100, .1, .1,
                           function(v) dnorm(v, log=T))
