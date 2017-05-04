## Based on example_logspec_partial
##setwd("~/research/gcp/hierarchical-estimation/logspec-bayes")

source("../example/logspec-data.R")

library(rstan)

stan.code <- "
data {
    int<lower=1> I; // number of observation
    int<lower=1> N; // number of ADM+
    int<lower=1> M; // number of ADM-
    int<lower=1> K; // number of coefficients (not including dropped)
    int<lower=1> L; // number of predictors, including intercept

    int<lower=1, upper=N> nn[I]; // ADM+ region for each observation
    int<lower=1, upper=M> mm[I]; // ADM- region for each observation

    vector[I] dmxx[K]; // demeaned predictors by observation
    matrix[I, L] zz[K]; // covariates per predictor
    vector[I] dmyy; // demeaned outcomes
}
parameters {
    real beta[K];
    vector[L] gamma[K]; // surface parameters

    vector<lower=0>[M] sigma;
}
model {
    { // Local scope for obs_mean
        vector[I] obs_mean;

        obs_mean = beta[1] * dmxx[1] .* exp(zz[1] * gamma[1]);
        for (kk in 2:K)
            obs_mean = obs_mean + beta[kk] * dmxx[kk] .* exp(zz[kk] * gamma[kk]);

        dmyy ~ normal(obs_mean, sigma[mm]);
    }
}"

stan.data$dmyy <- regional.demean(stan.data$yy, stan.data$nn)
stan.data$dmxx <- stan.data$xx
for (kk in 1:stan.data$K)
    stan.data$dmxx[kk, ] <- regional.demean(stan.data$dmxx[kk, ], stan.data$nn)

## Apply multiple times
iters <- 1000
portion.block <- 1/100
per.adm1s <- max(ceiling(portion.block * stan.data$M), stan.data$K * stan.data$L + 2)
portion.obs <- 1 - (per.adm1s - ceiling(portion.block * stan.data$M)) / per.adm1s
fit <- NULL

betas <- matrix(NA, 0, stan.data$K)
gammas <- matrix(NA, 0, stan.data$K * stan.data$L)

for (ii in 1:iters) {
    my.adm1s <- sample(1:stan.data$M, per.adm1s)
    included <- which(stan.data$mm %in% my.adm1s)
    if (portion.obs < 1)
        included <- sample(included, ceiling(length(included) * portion.obs))

    ## Construct bootstrapped dataset
    my.stan.data <- list(I = length(included),
                         N = length(unique(stan.data$nn[included])),
                         M = length(my.adm1s), K = stan.data$K, L = stan.data$L,
                         nn = as.numeric(factor(stan.data$nn[included])),
                         mm = as.numeric(factor(stan.data$mm[included])),
                         dmxx = stan.data$dmxx[, included],
                         zz = stan.data$zz[, included, ],
                         dmyy = stan.data$dmyy[included])

    ## Fit the model
    if (is.null(fit)) {
        fit <- stan(model_code = stan.code, data = my.stan.data,
                    iter = 1000, chains = 4)
    } else {
        fit <- stan(fit = fit, data = my.stan.data,
                    iter = 1000, chains = 4)
    }

    la <- extract(fit, permute=T)
    betas <- rbind(betas, matrix(colMeans(la$beta), 1, stan.data$K))
    gammas <- rbind(gammas, matrix(colMeans(la$gamma), 1, stan.data$K * stan.data$L))
}

