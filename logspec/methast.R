## Single Metropolis-Hastings with automatic tuning
## A 44% acceptance rate is optimal for 1-D sampling
## Want (a^44) * (b^56) = 1; Say a = 1.1, (1 / (1.1^44))^(1 / 56) = 0.9278487
methast <- function(iter, param0, sd0, calc.likeli) {
    params = matrix(NA, iter, length(param0))
    params[1, ] = param0

    last.likeli <- calc.likeli(param0)

    ## Look out for an even better solution
    best.likeli <- last.likeli
    best.index <- 1

    sd.product <- 1

    for (ii in 2:iter) {
        if (ii %% 100 == 0)
            print(ii)

        param.sample <- rnorm(length(param0), params[ii-1,], sd0 * sd.product)
        this.likeli <- calc.likeli(param.sample)

        prob <- exp(this.likeli - last.likeli)

        if (min(prob, 1) > runif(1)) {
            params[ii, ] <- param.sample
            last.likeli <- this.likeli
            sd.product <- sd.product * 1.1 # was too modest

            if (this.likeli > best.likeli) {
                best.likeli <- this.likeli
                best.index <- ii
            }
        } else {
            params[ii, ] <- params[ii-1, ]
            sd.product <- sd.product * 0.9278487 # was too bold
        }
    }

    list(params=params, best.index=best.index, best.likeli=best.likeli)
}

## Use Metropolis-Hastings with N seeds
repeated.methast <- function(seeds, iter, warmup, param0, sd0, likeli) {
    params <- matrix(NA, 0, length(param0))

    for (seed in 1:seeds) {
        print(paste("Seed", seed))
        methast.result <- methast(iter, param0, sd0, likeli)

        params <- rbind(params, as.matrix(methast.result$params[(warmup+1):iter,], iter-warmup, length(param0)))

        if (methast.result$best.index != 1) {
            param0 <- methast.result$params[methast.result$best.index,]
        }
    }

    list(params=params, best.param=param0)
}

parallel.single.methast <- function(prefix, seed, betas, gammas, betaerr, gammaerr, yy, xxs, zzs, adm1, adm2, weights=1, iter=600, warmup=100) {
    list2env(check.arguments(yy, xxs, zzs, adm1, adm2), parent.frame())
    list2env(demean.yxs(yy, xxs, adm2), parent.frame())

    methast.result <- methast.betagamma(K, L, dmxxs, dmyy, zzs, adm1, iter, betas, gammas, betaerr, gammaerr, weights)
    save(methast.result, file=paste0("MH-", prefix, seed, ".RData"))
}

parallel.combine.methast <- function(prefix, seeds, yy, xxs, zzs, adm1, adm2, weights=1, iter=600, warmup=100) {
    list2env(check.arguments(yy, xxs, zzs, adm1, adm2), parent.frame())
    list2env(demean.yxs(yy, xxs, adm2), parent.frame())

    betas <- matrix(NA, 0, K)
    gammas <- matrix(NA, 0, K*L)

    best.likeli <- -Inf
    beta0 <- NULL
    gamma0 <- NULL

    for (seed in 1:seeds) {
        if (!file.exists(paste0("MH-", prefix, seed, ".RData")))
            next

        load(paste0("MH-", prefix, seed, ".RData"))

        betas <- rbind(betas, methast.result$betas[(warmup+1):iter,])
        gammas <- rbind(gammas, methast.result$gammas[(warmup+1):iter,])

        this.likeli <- calc.likeli.nosigma(K, L, dmxxs, dmyy, zzs, adm1, methast.result$betas[methast.result$best.index,], matrix(methast.result$gammas[methast.result$best.index,], K, L), weights)
        if (this.likeli > best.likeli) {
            beta0 <- methast.result$betas[methast.result$best.index,]
            gamma0 <- matrix(methast.result$gammas[methast.result$best.index,], K, L)
        }
    }

    list(betas=betas, gammas=gammas, best.beta=beta0, best.gamma=gamma0)
}
