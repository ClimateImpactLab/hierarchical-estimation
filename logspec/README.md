# Using `estimate.logspec`

Please see `test/logspec/test_search.R` for an example.

## Terminology

Let there be `N` observations, each described by `K` weather
predictors (e.g., 11 bins, excluding a dropped bin).  The marginal
effect of each predictor varies according to up to `L` covariates
(e.g., 3 for mean temperature, log GDP p.c., and log population
density).

Each observation is associated with an ADM2 subregion of an ADM1
region.  There are `M` ADM1 regions.

## Arguments
 - `yy`: a vector of size `N`, containing the response for each observation
 - `xxs`: a `data.frame` or matrix of size `N x K`, containing predictors.
 - `zzs`: a `data.frame` or matrix of size `M x L`, containing covariates.  The rows should correspond to the values in `adm1` below (so, row 1 should give the covariates for administrative region 1).
 - `kls`: a matrix of size `K x L` of true and false.  Row `k`, column `l` of `kls` should be true if predictor `k` is informed by covariate `l`.
 - `adm1`: a vector of size `N`, containing values from 1 to `M`, attributing each observation to an ADM1 region.
 - `adm2`: a vector of size `N`, containing values from 1 to the total number of ADM2 subregions across all ADM1 regions.

Note that `adm1` can be generated from a collection of state names in
`data$state` using `as.numeric(factor(data$state))`, and `adm2` can be
generated from a collection of stat and county names using
`as.numeric(factor(paste(data$state, data$county)))`.

## Calling and result

First, you need to load the library:
```
source(".../logspec/search.R", chdir=T)
```

Call the estimate function as follows:
```
result <- estimate.logspec(yy, xxs, zzs, adm1, adm2)
```

The function will report each iteration it performs, and the
corresponding log likelihood.  After the iterations converge, the
numerical Hessian will be calculated to estimate the standard errors.

The result, `result`, is a list containing four sets of values:
 - `betas`: a vector of length `K`, containing the predictor coefficients.
 - `gammas`: a vector of size `sum(kls)`, containing the covariate coefficients.  The order is by predictor.
 - `sigma`: a vector of heteroskedastic error standard deviations for each ADM2 region.
 - `ses.betas`: A vector of length `K`, containing the standard error for each predictor coefficient.
 - `ses.gammas`: A matrix of size `K x L`, containing the standard error for each covariate coefficient.

## Table-based API

For datasets that are stored as tables of observations, the
table-based API wraps the main functionality of the system for more
convenience.  Two functions are provided, `ta.estimate.logspec` to
call `estimate.logspec` and `ta.estimate.vcv` to call `estimate.vcv`.
These functions replace the parameters `yy, xxs, zzs, kls, adm1, adm2`
with more table-friendly `df, outname, adm1name, adm2name, prednames,
covarnames`, as defined below:

 - `df`: A data.frame with all of the relevant information.
 - `outname`: A string for the name of the column describing the outcome variable.
 - `adm1name` and `adm2name`: String for the names of the columns defining the ADM1 and ADM2 regions.
 - `prednames`: The names of the columns for the predictors.  Each predictor name is repeated as many times as the covariates it has, plus 1 for an intercept covariate.
 - `covarnames`: The names of the columns for the covariates.  This should have as many entries as the `prednames` vector, where the corresponding entry for a given predictor should be one of its covariates.  Each predictor should have an entry for the covariate '1'.

For example, the 5-bin example is called as follows:
```
source(".../logspec/tableapi.R", chdir=T)
ta.estimate.logspec(df, 'rate', 'adm1', 'adm2',
                    c('bin1', 'bin1', 'bin1', 'bin2', 'bin2', 'bin2',
                      'bin4', 'bin4', 'bin4', 'bin5', 'bin5', 'bin5'),
                    c('1', 'meant', 'log_gdppc', '1', 'meant', 'log_gdppc',
                      '1', 'meant', 'log_gdppc', '1', 'meant', 'log_gdppc'))
```
