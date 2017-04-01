# Using `estimate.logspec`

Please see `example.R` for an example.

## Terminology

Let there be `N` observations, each described by `K` weather
predictors (e.g., 11 bins, excluding a dropped bin).  The marginal
effect of each predictor varies according to `L` covariates (e.g., 3
for mean temperature, log GDP p.c., and log population density).

Each observation is associated with an ADM2 subregion of an ADM1
region.  There are `M` ADM1 regions.

## Arguments
 - `yy`: a vector of size `N`, containing the response for each observation
 - `xxs`: a `data.frame` or matrix of size `N x K`, containing predictors.
 - `zzs`: a `data.frame` or matrix of size `M x KL`, containing covariates.  The rows should correspond to the values in `adm1` below (so, row 1 should give the covariates for administrative region 1).  The columns should be organized so that the covariates for a given predictor are grouped together (so, columns 1 - 4 might be all used to describe the marginal response of days in one temperature bin; columns 5 - 9 would be another bin; and so on).
 - `adm1`: a vector of size `N`, containing values from 1 to `M`, attributing each observation to an ADM1 region.
 - `adm2`: a vector of size `N`, containing values from 1 to the total number of ADM2 subregions across all ADM1 regions.

Note that `adm1` can be generated from a collection of state names in
`data$state` using `as.numeric(factor(data$state))`, and `adm2` can be
generated from a collection of stat and county names using
`as.numeric(factor(paste(data$state, data$county)))`.

Given a dataset with covariates for every observation, the `zzs`
matrix can be formed with `data[!duplicated(adm1), COLUMNS]`.

## Calling and result

Call the estimate function as follows:
```
result <- estimate.logspec(yy, xxs, zzs, adm1, adm2)
```

The function will report each iteration it performs, and the
corresponding log likelihood.  After the iterations converge, the
numerical Hessian will be calculated to estimate the standard errors.

The result, `result`, is a list containing four sets of values:
 - `betas`: a vector of length `K`, containing the predictor coefficients
 - `gammas`: a matrix of size `K x L`, containing the covariate coefficients
 - `ses.betas`: A vector of length `K`, containing the standard error for each predictor coefficient
 - `ses.gammas`: A matrix of size `K x L`, containing the standard error for each covariate coefficient

