##setwd("~/research/gcp/hierarchical-estimation/logspec")

library(rstan)

df <- read.csv("binned.csv")

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

## Prepare the data
df$log_gdppc <- log(df$gdppc)
df$adm1 <- factor(df$adm1)
df$adm2 <- factor(df$adm2)

stan.data <- list(I = nrow(df), N = length(unique(df$adm2)),
                  M = length(unique(df$adm1)), K = 4, L = 2,
                  nn = as.numeric(df$adm2), mm = as.numeric(df$adm1),
                  xx = t(df[, c('bin1', 'bin2', 'bin4', 'bin5')]),
                  zz = aperm(array(unlist(list(df[, c('meant', 'log_gdppc')],
                                               df[, c('meant', 'log_gdppc')],
                                               df[, c('meant', 'log_gdppc')],
                                               df[, c('meant', 'log_gdppc')])),
                                   dim = c(nrow(df), 2, 4)), c(3, 1, 2)),
                  yy = df$rate)

## Fit the model
fit <- stan(model_code = stan.code, data = stan.data,
            iter = 1000, chains = 4)
