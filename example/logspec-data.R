df <- read.csv("binned.csv")

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
