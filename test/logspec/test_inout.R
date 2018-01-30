library(RUnit)

## Load the library
source("logspec/tableapi.R", chdir=T)

tstar <- 22

test.ta.estimate.logspec <- function() {
    ## Load the data
    df <- read.csv("example/true-binned.csv")
    df$log_gdppc <- log(df$gdppc)

    result <- ta.estimate.logspec(df, 'rate', 'adm1',
                                  c('bin1', 'bin1', 'bin1', 'bin2', 'bin2', 'bin2',
                                    'bin4', 'bin4', 'bin4', 'bin5', 'bin5', 'bin5'),
                                  c('1', 'meant', 'log_gdppc', '1', 'meant', 'log_gdppc',
                                    '1', 'meant', 'log_gdppc', '1', 'meant', 'log_gdppc'), 'adm2')

    vcvfit <- ta.estimate.vcv(result$betas, result$gammas, result$sigma,
                              df, 'rate', 'adm1',
                              c('bin1', 'bin1', 'bin1', 'bin2', 'bin2', 'bin2',
                                'bin4', 'bin4', 'bin4', 'bin5', 'bin5', 'bin5'),
                              c('1', 'meant', 'log_gdppc', '1', 'meant', 'log_gdppc',
                                '1', 'meant', 'log_gdppc', '1', 'meant', 'log_gdppc'), 'adm2')

    ta.output.csvv(df, 'rate', 'adm1',
                   c('bin1', 'bin1', 'bin1', 'bin2', 'bin2', 'bin2',
                     'bin4', 'bin4', 'bin4', 'bin5', 'bin5', 'bin5'),
                   c('1', 'meant', 'log_gdppc', '1', 'meant', 'log_gdppc',
                     '1', 'meant', 'log_gdppc', '1', 'meant', 'log_gdppc'),
                   'adm2', result, vcvfit,
                   "This is a test.",
                   "TEST-20180129", c(),
                   "This is a longer description.",
                   list(rate="Rate of something", bin1="Temperature bin 1", bin2="Temperature bin 2",
                        bin4="Temperature bin 4", bin5="Temperature bin 5", meant="Mean temperature",
                        log_gdppc="Log GDP per capita"),
                   list(rate="widgets / s", bin1="C", bin2="C", bin4="C", bin5="C", meant="C", log_gdppc="log USD"),
                   "output.csvv")
}
