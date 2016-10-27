"
Stage 2 mortality regressions

:author: My Name
:contact: `myname@gmail.com <mailto:myname@gmail.com>`_
:modified: 2016/08/19
:team: Mortality
:lead: Amir Jina

Mortality Stage 2
* Stage 1 : 5C bins, 13-18 reference, population weighted temperature, 'unweighted' (adm2=1) regression
* Stage 2 : Regress Stage 1 coefficients on:
  - average # days in bin (population weighted, averaged over sample period)
  - average GDP per capita, log
  - average population density, log
  - average share of population age < 5
  - average share of population age > 65

options for excluding outliers, making covariates cubic, and (un)weighting

Note: This is an EXAMPLE header and was created to look real-ISH! DO NOT 
INTERPRET THIS INFORMATION AS ACCURATE!!!


Include files
~~~~~~~~~~~~~
surface.R
    :location: ~/Dropbox/ChicagoRAs_Reanalysis/Repos/impact-calculations/interpolate/surface.R
    :description: build interpolation surface
    :version: GCP.Mortality.surface.2016-08-19


Input variables
~~~~~~~~~~~~~~~

intergroup_tp3_08_31_2016_manual.csv
    :location: ~/Dropbox/conflict/Conflict_reanalysis/output/Final/stage1
    :description: Intergroup conflict regressed on T + P^3
    :version: ``GCP.Conflict.intergroup_tp3_manual.2016-08-31``


Output variables
~~~~~~~~~~~~~~~~

intergroup_fit_objects.rda
    :location: ~/Dropbox/GCP/CONFLICT/results/stage2
    :description: intergroup conflict fit
    :version: ``GCP.Conflict.intergroup_fit_objects.2016-08-31``

"

# PRELIMINARIES ============================================

# Location of this file
setwd("~/Dropbox/conflict/Conflict_reanalysis/code")

# Load James's estimation code, let Stan use multicores
#source("~/research/gcp/impact-calculations/interpolate/surface.R")
source("~/Dropbox/ChicagoRAs_Reanalysis/Repos/impact-calculations/interpolate/surface.R")
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

# Date for output file versioning
datesuff <- format(Sys.Date(), "%Y_%m_%d")
# Input/output location and
dir.input <- "~/Dropbox/conflict/Conflict_reanalysis/output/Final/stage1"
dir.output <- "~/Dropbox/GCP/CONFLICT/results/stage2"
dir.tables <- "~/Dropbox/GCP/CONFLICT/results/tables"
coeffile <- paste0(dir.input, "/intergroup_tp3_08_31_2016_manual.csv")

# ESTIMATION ===============================================

# Unweighted, cubic 2nd stage covariates -------------------

# Set specification options and name for output files
do.nooutliers <- T
do.unweighted <- T
do.cubicpreds <- T
spec <- "intergroup_tp3_unw_cub"

# 1st Stage climate variables
climvars <- c('temp', 'precip1', 'precip2', 'precip3')
K <- length(climvars)
coefvars <- paste0('beta_', climvars)
serrvars <- paste0('se_', climvars)

# 2nd Stage covariates
predvars <- c('avgT', 'avgprecip', 'logpopop', 'loggdppc')
L.nocons <- length(predvars)
if (do.cubicpreds) {
  L <- 1 + L.nocons*3
  allpreds <- c('const', predvars, paste0(predvars, 2), paste0(predvars, 3))
} else {
  L <- 1 + L.nocons
  allpreds <- c('const', predvars)
}

# Create a BayesObservations object to hold the data
surface <- SurfaceObservations(K=K, L=L)

# Read inputs
inputs <- read.csv(coeffile)
inputs$loggdppc <- log(inputs$gdppc)
inputs$logpopop <- log(inputs$popop)
# Drop outliers
if (do.nooutliers) {
  # Based on temperature coefficients
  zscore <- (inputs$beta_temp - mean(inputs$beta_temp)) / sd(inputs$beta_temp)
  inputs <- subset(inputs, abs(zscore) < 4)
  zscore <- (inputs$beta_temp - mean(inputs$beta_temp)) / sd(inputs$beta_temp)
  inputs <- subset(inputs, abs(zscore) < 4)
  zscore <- (inputs$beta_temp - mean(inputs$beta_temp)) / sd(inputs$beta_temp)
  inputs <- subset(inputs, abs(zscore) < 4)
  
  # Based on 1st order precipitation coefficients
  ss <- subset(inputs, !is.na(inputs$beta_precip1))
  zscore <- (ss$beta_precip1 - mean(ss$beta_precip1)) / sd(ss$beta_precip1)
  inputs <- subset(inputs, abs(zscore) < 4)
}

# EXAMPLE of reading in stage 1 VCV matrix
#vcvs <- read.csv("stage1_output_interpersonal_property_vcvmatrices.csv")

# Add each row of stage 1 output as a new observation
for (jj in 1:nrow(inputs)) {
  # Get coefficients and standard errors from input file
  betas <- inputs[jj, coefvars]
  betas[is.na(betas)] <- 0
  names(betas) <- climvars
  serrs <- inputs[jj, serrvars]
  serrs[is.na(serrs)] <- Inf

  # TEMPORARY - Drop tiny betas/SEs to avoid positive definite error
  # betas[serrs < 1e-4] <- 0
  # serrs[serrs < 1e-4] <- Inf
  
  # If VCVs are not available, just diagonalize SEs
  if (do.unweighted)
      vcv <- diag(rep(1, length(serrs)))
  else
      vcv <- diag(as.numeric(serrs)^2)
  names(vcv) <- climvars

  # If VCVs are available, read them
  # rr <- inputs[jj, 'region_code']
  # vcv <- vcvs[vcvs$region_code==rr, climvars]

  # 
  predses <- matrix(0, 0, L)
  for (kk in 1:K) {
      if (do.cubicpreds) {
          row <- cbind(data.frame(const=1), inputs[jj, predvars], inputs[jj, predvars]^2, inputs[jj, predvars]^3)
      } else
          row <- cbind(data.frame(const=1), inputs[jj, predvars])
      predses <- rbind(predses, row)
  }
  names(predses) <- allpreds

  # Make sure that all SEs haven't been dropped
  if (min(serrs == Inf) == 0 && max(is.na(predses)) == 0) {
    surface <- addObs(surface, betas, vcv, predses)
  }

}

# Seemingly Unrelated Regression (SUR) .....................

fit.sur <- estimate.semur(surface)
# summary(fit.sur)

## Output bin surface parameters
## dependencies in a single file?
surface.write(surface, fit.sur, paste0(dir.output, "/", spec, "_semur_auto_", datesuff, ".csvv"), "Conflict group temp stage 2 results", "CONFLICT-INTERGROUP-TP3", "intergroup_tp3.csv", climvars, allpreds)


# Bayesian (BAY) ...........................................

fit.bay <- estimate.bayes(surface)
# print(fit.bay)
#la <- extract(fit.bay, permute=T)


## Output bin surface parameters
## dependencies in a single file?
surface.write(surface, fit.bay, paste0(dir.output, "/", spec, "_bayes_auto_", datesuff, ".csvv"), "Conflict group temp stage 2 results", "CONFLICT-INTERGROUP-TP3", "intergroup_tp3.csv", climvars, allpreds)


fit.sur.cub.unw <- fit.sur
fit.bay.cub.unw <- fit.bay


# Variance weighted, linear 2nd stage covariates -----------

do.nooutliers <- T
do.unweighted <- F
do.cubicpreds <- F
spec <- "intergroup_tp3_var_lin"

## Specification of 1st Stage model
K <- 4

## Specification of 2nd Stage model
if (do.cubicpreds) {
  L <- 1 + 4*3
} else
  L <- 5

## Create a BayesObservations object to hold the data
surface <- SurfaceObservations(K=K, L=L)

## Read inputs
climvars <- c('temp', 'precip1', 'precip2', 'precip3')
coefvars <- paste0('beta_', climvars)
serrvars <- paste0('se_', climvars)
predvars <- c('avgT', 'avgprecip', 'logpopop', 'loggdppc')
if (do.cubicpreds) {
  allpreds <- c('const', predvars, paste0(predvars, 2), paste0(predvars, 3))
} else
  allpreds <- c('const', predvars)

##inputs <- read.csv("stage1_output_interpersonal_property_coefficients.csv")

inputs <- read.csv(coeffile)
inputs$loggdppc <- log(inputs$gdppc)
inputs$logpopop <- log(inputs$popop)
if (do.nooutliers) {
  zscore <- (inputs$beta_temp - mean(inputs$beta_temp)) / sd(inputs$beta_temp)
  inputs <- subset(inputs, abs(zscore) < 4)
  zscore <- (inputs$beta_temp - mean(inputs$beta_temp)) / sd(inputs$beta_temp)
  inputs <- subset(inputs, abs(zscore) < 4)
  zscore <- (inputs$beta_temp - mean(inputs$beta_temp)) / sd(inputs$beta_temp)
  inputs <- subset(inputs, abs(zscore) < 4)
  
  ss <- subset(inputs, !is.na(inputs$beta_precip1))
  zscore <- (ss$beta_precip1 - mean(ss$beta_precip1)) / sd(ss$beta_precip1)
  inputs <- subset(inputs, abs(zscore) < 4)
}

#vcvs <- read.csv("stage1_output_interpersonal_property_vcvmatrices.csv")

for (jj in 1:nrow(inputs)) {
  betas <- inputs[jj, coefvars]
  betas[is.na(betas)] <- 0
  names(betas) <- climvars
  
  ## If VCVs are not available, just diagonalize SEs
  serrs <- inputs[jj, serrvars]
  serrs[is.na(serrs)] <- Inf
  
  # TEMPORARY - Drop tiny betas/SEs to avoid positive definite error
  #betas[serrs < 1e-4] <- 0
  #serrs[serrs < 1e-4] <- Inf
  
  if (do.unweighted)
    vcv <- diag(rep(1, length(serrs)))
  else
    vcv <- diag(as.numeric(serrs)^2)
  
  names(vcv) <- climvars
  
  ## If VCVs are available, read them
  #rr <- inputs[jj, 'region_code']
  #vcv <- vcvs[vcvs$region_code==rr, climvars]
  
  predses <- matrix(0, 0, L)
  for (kk in 1:length(predvars)) {
    if (do.cubicpreds) {
      row <- cbind(data.frame(const=1), inputs[jj, predvars], inputs[jj, predvars]^2, inputs[jj, predvars]^3)
    } else
      row <- cbind(data.frame(const=1), inputs[jj, predvars])
    predses <- rbind(predses, row)
  }
  names(predses) <- allpreds
  
  # Make sure that all SEs haven't been dropped
  if (min(serrs == Inf) == 0 && max(is.na(predses)) == 0) {
    surface <- addObs(surface, betas, vcv, predses)
  }
  
}

## Test SUR
fit.sur <- estimate.semur(surface)
# summary(fit.sur)

## Output bin surface parameters
# dependencies in a single file?
surface.write(surface, fit.sur, paste0(dir.output, "/", spec, "_semur_auto_", datesuff, ".csvv"), "Conflict group temp stage 2 results", "CONFLICT-INTERGROUP-TP3", "intergroup_tp3.csv", climvars, allpreds)


## Test Bayes
fit.bay <- estimate.bayes(surface)
# print(fit.bay)

#la <- extract(fit.bay, permute=T)

## Output bin surface parameters
## dependencies in a single file?
surface.write(surface, fit.bay, paste0(dir.output, "/", spec, "_bayes_auto_", datesuff, ".csvv"), "Conflict group temp stage 2 results", "CONFLICT-INTERGROUP-TP3", "intergroup_tp3.csv", climvars, allpreds)

fit.sur.lin.var <- fit.sur
fit.bay.lin.var <- fit.bay


# Variance weighted, cubic 2nd stage covariates ------------

do.nooutliers <- T
do.unweighted <- F
do.cubicpreds <- T
spec <- "intergroup_tp3_var_cub"

## Specification of 1st Stage model
K <- 4

## Specification of 2nd Stage model
if (do.cubicpreds) {
  L <- 1 + 4*3
} else
  L <- 5

## Create a BayesObservations object to hold the data
surface <- SurfaceObservations(K=K, L=L)

## Read inputs
climvars <- c('temp', 'precip1', 'precip2', 'precip3')
coefvars <- paste0('beta_', climvars)
serrvars <- paste0('se_', climvars)
predvars <- c('avgT', 'avgprecip', 'logpopop', 'loggdppc')
if (do.cubicpreds) {
  allpreds <- c('const', predvars, paste0(predvars, 2), paste0(predvars, 3))
} else
  allpreds <- c('const', predvars)

##inputs <- read.csv("stage1_output_interpersonal_property_coefficients.csv")

inputs <- read.csv(coeffile)
inputs$loggdppc <- log(inputs$gdppc)
inputs$logpopop <- log(inputs$popop)
if (do.nooutliers) {
  zscore <- (inputs$beta_temp - mean(inputs$beta_temp)) / sd(inputs$beta_temp)
  inputs <- subset(inputs, abs(zscore) < 4)
  zscore <- (inputs$beta_temp - mean(inputs$beta_temp)) / sd(inputs$beta_temp)
  inputs <- subset(inputs, abs(zscore) < 4)
  zscore <- (inputs$beta_temp - mean(inputs$beta_temp)) / sd(inputs$beta_temp)
  inputs <- subset(inputs, abs(zscore) < 4)
  
  ss <- subset(inputs, !is.na(inputs$beta_precip1))
  zscore <- (ss$beta_precip1 - mean(ss$beta_precip1)) / sd(ss$beta_precip1)
  inputs <- subset(inputs, abs(zscore) < 4)
}

#vcvs <- read.csv("stage1_output_interpersonal_property_vcvmatrices.csv")


for (jj in 1:nrow(inputs)) {
  betas <- inputs[jj, coefvars]
  betas[is.na(betas)] <- 0
  names(betas) <- climvars
  
  ## If VCVs are not available, just diagonalize SEs
  serrs <- inputs[jj, serrvars]
  serrs[is.na(serrs)] <- Inf
  
  # TEMPORARY - Drop tiny betas/SEs to avoid positive definite error
  #betas[serrs < 1e-4] <- 0
  #serrs[serrs < 1e-4] <- Inf
  
  if (do.unweighted)
    vcv <- diag(rep(1, length(serrs)))
  else
    vcv <- diag(as.numeric(serrs)^2)
  
  names(vcv) <- climvars
  
  ## If VCVs are available, read them
  #rr <- inputs[jj, 'region_code']
  #vcv <- vcvs[vcvs$region_code==rr, climvars]
  
  predses <- matrix(0, 0, L)
  for (kk in 1:length(predvars)) {
    if (do.cubicpreds) {
      row <- cbind(data.frame(const=1), inputs[jj, predvars], inputs[jj, predvars]^2, inputs[jj, predvars]^3)
    } else
      row <- cbind(data.frame(const=1), inputs[jj, predvars])
    predses <- rbind(predses, row)
  }
  names(predses) <- allpreds
  
  # Make sure that all SEs haven't been dropped
  if (min(serrs == Inf) == 0 && max(is.na(predses)) == 0) {
    surface <- addObs(surface, betas, vcv, predses)
  }
  
}

## Test SUR

fit.sur <- estimate.semur(surface)
# summary(fit.sur)

## Output bin surface parameters
## dependencies in a single file?
surface.write(surface, fit.sur, paste0(dir.output, "/", spec, "_semur_auto_", datesuff, ".csvv"), "Conflict group temp stage 2 results", "CONFLICT-INTERGROUP-TP3", "intergroup_tp3.csv", climvars, allpreds)



## Test Bayes
fit.bay <- estimate.bayes(surface)
# print(fit.bay)
#la <- extract(fit.bay, permute=T)


## Output bin surface parameters
## dependencies in a single file?
surface.write(surface, fit.bay, paste0(dir.output, "/", spec, "_bayes_auto_", datesuff, ".csvv"), "Conflict group temp stage 2 results", "CONFLICT-INTERGROUP-TP3", "intergroup_tp3.csv", climvars, allpreds)

fit.sur.cub.var <- fit.sur
fit.bay.cub.var <- fit.bay


# Unweighted, linear 2nd stage covariates ------------------

do.nooutliers <- T
do.unweighted <- T
do.cubicpreds <- F
spec <- "intergroup_tp3_unw_lin"

## Specification of 1st Stage model
K <- 4

## Specification of 2nd Stage model
if (do.cubicpreds) {
  L <- 1 + 4*3
} else
  L <- 5

## Create a BayesObservations object to hold the data
surface <- SurfaceObservations(K=K, L=L)

## Read inputs
climvars <- c('temp', 'precip1', 'precip2', 'precip3')
coefvars <- paste0('beta_', climvars)
serrvars <- paste0('se_', climvars)
predvars <- c('avgT', 'avgprecip', 'logpopop', 'loggdppc')
if (do.cubicpreds) {
  allpreds <- c('const', predvars, paste0(predvars, 2), paste0(predvars, 3))
} else
  allpreds <- c('const', predvars)

##inputs <- read.csv("stage1_output_interpersonal_property_coefficients.csv")

inputs <- read.csv(coeffile)
inputs$loggdppc <- log(inputs$gdppc)
inputs$logpopop <- log(inputs$popop)
if (do.nooutliers) {
  zscore <- (inputs$beta_temp - mean(inputs$beta_temp)) / sd(inputs$beta_temp)
  inputs <- subset(inputs, abs(zscore) < 4)
  zscore <- (inputs$beta_temp - mean(inputs$beta_temp)) / sd(inputs$beta_temp)
  inputs <- subset(inputs, abs(zscore) < 4)
  zscore <- (inputs$beta_temp - mean(inputs$beta_temp)) / sd(inputs$beta_temp)
  inputs <- subset(inputs, abs(zscore) < 4)
  
  ss <- subset(inputs, !is.na(inputs$beta_precip1))
  zscore <- (ss$beta_precip1 - mean(ss$beta_precip1)) / sd(ss$beta_precip1)
  inputs <- subset(inputs, abs(zscore) < 4)
}

#vcvs <- read.csv("stage1_output_interpersonal_property_vcvmatrices.csv")


for (jj in 1:nrow(inputs)) {
  betas <- inputs[jj, coefvars]
  betas[is.na(betas)] <- 0
  names(betas) <- climvars
  
  ## If VCVs are not available, just diagonalize SEs
  serrs <- inputs[jj, serrvars]
  serrs[is.na(serrs)] <- Inf
  
  # TEMPORARY - Drop tiny betas/SEs to avoid positive definite error
  #betas[serrs < 1e-4] <- 0
  #serrs[serrs < 1e-4] <- Inf
  
  if (do.unweighted)
    vcv <- diag(rep(1, length(serrs)))
  else
    vcv <- diag(as.numeric(serrs)^2)
  
  names(vcv) <- climvars
  
  ## If VCVs are available, read them
  #rr <- inputs[jj, 'region_code']
  #vcv <- vcvs[vcvs$region_code==rr, climvars]
  
  predses <- matrix(0, 0, L)
  for (kk in 1:length(predvars)) {
    if (do.cubicpreds) {
      row <- cbind(data.frame(const=1), inputs[jj, predvars], inputs[jj, predvars]^2, inputs[jj, predvars]^3)
    } else
      row <- cbind(data.frame(const=1), inputs[jj, predvars])
    predses <- rbind(predses, row)
  }
  names(predses) <- allpreds
  
  # Make sure that all SEs haven't been dropped
  if (min(serrs == Inf) == 0 && max(is.na(predses)) == 0) {
    surface <- addObs(surface, betas, vcv, predses)
  }
  
}

## Test SUR

fit.sur <- estimate.semur(surface)
# summary(fit.sur)

## Output bin surface parameters
## dependencies in a single file?
surface.write(surface, fit.sur, paste0(dir.output, "/", spec, "_semur_auto_", datesuff, ".csvv"), "Conflict group temp stage 2 results", "CONFLICT-INTERGROUP-TP3", "intergroup_tp3.csv", climvars, allpreds)



## Test Bayes
fit.bay <- estimate.bayes(surface)
# print(fit.bay)
#la <- extract(fit.bay, permute=T)


## Output bin surface parameters
## dependencies in a single file?
surface.write(surface, fit.bay, paste0(dir.output, "/", spec, "_bayes_auto_", datesuff, ".csvv"), "Conflict group temp stage 2 results", "CONFLICT-INTERGROUP-TP3", "intergroup_tp3.csv", climvars, allpreds)

fit.sur.lin.unw <- fit.sur
fit.bay.lin.unw <- fit.bay

save(fit.sur.lin.unw, fit.sur.lin.var, fit.sur.cub.unw, fit.sur.cub.var, fit.bay.lin.unw, fit.bay.lin.var, fit.bay.cub.unw, fit.bay.cub.var,
     file = paste0(dir.output, "/intergroup_fit_objects.rda"))


# FORMATTED OUTPUT =========================================

# texred uses a different `extract` command to
library(texreg)

tr.sur.lin.var <- extract(fit.sur.lin.var)
tr.sur.cub.var <- extract(fit.sur.cub.var)
tr.sur.lin.unw <- extract(fit.sur.lin.unw)
tr.sur.cub.unw <- extract(fit.sur.cub.unw)
slu <- summary(fit.bay.lin.unw, pars="gamma")$summary
scu <- summary(fit.bay.cub.unw, pars="gamma")$summary
slv <- summary(fit.bay.lin.var, pars="gamma")$summary
scv <- summary(fit.bay.cub.var, pars="gamma")$summary
nl <- nrow(slv) / K
nc <- nrow(scu) / K
for (eq in 1:K) {
  # Convert Bayesian results to texreg format
  prednames <- names(fit.sur.cub.unw$eq[[eq]]$coefficients)
  ncoef <- length(prednames)
  
  s <- slu
  n <- nl
  coef <- vector(length=ncoef)
  ix1 <- 1+n*(eq-1)
  ix2 <- n*eq
  coef[1:n] <- s[ix1:ix2,"mean"]
  coef[-(1:n)] <- NA
  se <- vector(length=ncoef)
  se[1:n] <- s[ix1:ix2,"se_mean"]
  se[-(1:n)] <- NA
  tr.bay.lin.unw <- createTexreg(coef.names = prednames, coef = coef, se = se)

  s <- slv
  n <- nl
  coef <- vector(length=ncoef)
  ix1 <- 1+n*(eq-1)
  ix2 <- n*eq
  coef[1:n] <- s[ix1:ix2,"mean"]
  coef[-(1:n)] <- NA
  se <- vector(length=ncoef)
  se[1:n] <- s[ix1:ix2,"se_mean"]
  se[-(1:n)] <- NA
  tr.bay.lin.var <- createTexreg(coef.names = prednames, coef = coef, se = se)

  s <- scu
  n <- nc
  coef <- vector(length=ncoef)
  ix1 <- 1+n*(eq-1)
  ix2 <- n*eq
  coef[1:n] <- s[ix1:ix2,"mean"]
  coef[-(1:n)] <- NA
  se <- vector(length=ncoef)
  se[1:n] <- s[ix1:ix2,"se_mean"]
  se[-(1:n)] <- NA
  tr.bay.cub.unw <- createTexreg(coef.names = prednames, coef = coef, se = se)

  s <- scv
  n <- nc
  coef <- vector(length=ncoef)
  ix1 <- 1+n*(eq-1)
  ix2 <- n*eq
  coef[1:n] <- s[ix1:ix2,"mean"]
  coef[-(1:n)] <- NA
  se <- vector(length=ncoef)
  se[1:n] <- s[ix1:ix2,"se_mean"]
  se[-(1:n)] <- NA
  tr.bay.cub.var <- createTexreg(coef.names = prednames, coef = coef, se = se)

  caption = paste("Estimated relationship between",climvars[eq],"1st stage coefficients and 2nd stage covariates")
  covars <- c("Constant", "Avg. Temp.", "Avg. Precip.", "Pop. Density", "GDP Per Capita", "(Avg. Temp)$^2$", "(Avg. Precip.)$^2$", "(Pop. Density)$^2$", "(GDP Per Capita)$^2$", "(Avg. Temp)$^3$", "(Avg. Precip.)$^3$", "(Pop. Density)$^3$", "(GDP Per Capita)$^3$")
  models <- c("Linear (Unweighted)","Linear (Weighted)","Cubic (Unweighted)","Cubic (Weighted)","Linear (Unweighted)","Linear (Weighted)","Cubic (Unweighted)","Cubic (Weighted)")
  texreg(list(tr.sur.lin.unw[[eq]], tr.sur.lin.var[[eq]], tr.sur.cub.unw[[eq]], tr.sur.cub.var[[eq]], tr.bay.lin.unw, tr.bay.lin.var, tr.bay.cub.unw, tr.bay.cub.var), 
         file = paste0(dir.tables,"/intergroup_robustness_table_",climvars[eq],"_",datesuff,".tex"),
         custom.model.names = models,
         custom.coef.names = covars,
         digits = 5,
         caption = caption)
  # SUR only
  texreg(list(tr.sur.lin.unw[[eq]], tr.sur.lin.var[[eq]], tr.sur.cub.unw[[eq]], tr.sur.cub.var[[eq]]), 
         file = paste0(dir.tables,"/intergroup_robustness_table_sur_",climvars[eq],"_",datesuff,".tex"),
         custom.model.names = models[1:4],
         custom.coef.names = covars,
         digits = 5,
         caption = paste(caption,"(Seemingly Unrelated estimation)"))
  # Bayesian only
  texreg(list(tr.bay.lin.unw, tr.bay.lin.var, tr.bay.cub.unw, tr.bay.cub.var), 
         file = paste0(dir.tables,"/intergroup_robustness_table_bay_",climvars[eq],"_",datesuff,".tex"),
         custom.model.names = models[1:4],
         custom.coef.names = covars,
         digits = 5,
         caption = paste(caption,"(Bayesian estimation)"))
}

detach("package:texreg")

