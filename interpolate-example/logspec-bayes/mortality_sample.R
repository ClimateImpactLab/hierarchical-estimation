library(foreign)
df <- read.dta("mortality_foranalysis_withcov_12.dta")

adm2s <- sample(unique(df$adm2), 1000)
subdf <- subset(df, adm2 %in% adm2s)

subdf1 <- subdf[subdf$adm2 %in% adm2s[1:50],]
subdf2 <- subdf[subdf$adm2 %in% adm2s[51:1000],]

subdf <- rbind(subdf1, subdf2[sample(1:nrow(subdf2), 50000 - nrow(subdf1)),])
write.csv(subdf, "mortality_sub50000.csv", row.names=F)
