## Create a BayesObservations object to hold the data
source("surface.R")
surface <- SurfaceObservations(K=7, L=3)

data.dir <- "~/Dropbox/GCP/MORTALITY/tables/interaction/single_stage/new/ster/AH_New/adm1/Output"

betas <- read.csv(file.path(data.dir, "betas_adm1.csv"))
columns <- c('b_GMFD_term0_NS', 'b_GMFD_term1_NS', 'b_GMFD_term2_NS', 'b_GMFD_term3_NS', 'b_GMFD_term4_NS', 'b_GMFD_term5_NS', 'b_GMFD_term6_NS')

for (ii in 1:nrow(betas)) { # 90 worksish, 91 doesn'tish (tau=500)
    if (is.na(betas$gdppc[ii]))
        next
    ## Find associated VCV
    vcv.file <- Sys.glob(file.path(data.dir, paste0("cov_*_", ii, ".csv")))
    vcv <- as.matrix(read.csv(vcv.file, header=F))[1:7, 1:7]

    if (any(eigen(vcv)$values < 0)) {
        ses <- betas[ii, sub("b_", "se_", columns)]
        if (any(is.na(ses)))
            next
        vcv <- diag(ses)^2
    }

    predses <- cbind(rep(1, 7), rep(log(betas$gdppc[ii]), 7), rep(betas$Tmean_GMFD[ii], 7))

    surface <- addObs(surface, betas[ii, columns], vcv, predses)
}
