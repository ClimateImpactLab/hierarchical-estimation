##setwd("~/research/gcp/hierarchical-estimation/ranges")

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
transformed parameters {
    vector[I] obs_mean;

    obs_mean = beta[1] * dmxx[1] .* exp(zz[1] * gamma[1]);
    for (kk in 2:K)
        obs_mean = obs_mean + beta[kk] * dmxx[kk] .* exp(zz[kk] * gamma[kk]);
}
model {
    dmyy ~ normal(obs_mean, sigma[mm]);
}"

stan.data$dmyy <- regional.demean(stan.data$yy, stan.data$nn)
stan.data$dmxx <- stan.data$xx
for (kk in 1:stan.data$K)
    stan.data$dmxx[kk, ] <- regional.demean(stan.data$dmxx[kk, ], stan.data$nn)

library(RcppArmadillo)

betas <- c(8.1, 2.0, 2.0, 8.1)
gammas <- matrix(c(rep(-.005, 4), rep(-0.05, 4)), 4, 2)

## Determine the sigmas
dmxx.exp <- stan.data$dmxx
for (kk in 1:stan.data$K)
    dmxx.exp[kk, ] <- dmxx.exp[kk, ] * exp(as.matrix(stan.data$zz[, ((kk-1)*stan.data$L + 1):(kk*stan.data$L)]) %*% gammas[kk, ])

sigmas <- c()
for (jj in 1:stan.data$M) {
    included <- stan.data$mm == jj

    residuals <- t(dmxx.exp[, included]) %*% betas - stan.data$dmyy[included]

    sigmas <- c(sigmas, sd(residuals))
}

init <- function() {
    list(beta=betas, gamma=gammas, sigma=sigmas)
}

## Fit the model
fit <- stan(model_code = stan.code, data = stan.data,
            iter = 1000, chains = 4, init=init)
