##setwd("~/research/gcp/hierarchical-estimation/logspec")

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

    vector[I] xx[K]; // predictors by observation
    matrix[I, L] zz[K]; // covariates per predictor
    vector[I] yy; // outcomes
}
parameters {
    vector[N] alpha;
    real beta[K];
    vector[L] gamma[K]; // surface parameters

    vector<lower=0>[M] sigma;
}
transformed parameters {
    vector[I] obs_mean;

    obs_mean = alpha[nn];
    for (kk in 1:K)
        obs_mean = obs_mean + beta[kk] * xx[kk] .* exp(zz[kk] * gamma[kk]);
}
model {
    yy ~ normal(obs_mean, sigma[mm]);
}"

## Fit the model
fit <- stan(model_code = stan.code, data = stan.data,
            iter = 1000, chains = 4)
