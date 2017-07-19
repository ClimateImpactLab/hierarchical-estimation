library(ggplot2)

ggstandard <- function(allpoints, xx) {
    ggplot(data.frame(x=xx, y=allpoints[2,], ymin=allpoints[1,], ymax=allpoints[3,]), aes(x, y, ymin=ymin, ymax=ymax)) +
        geom_line() + geom_ribbon(alpha=.5) +
        xlab("Temperature") + ylab("Regular response") + theme_bw() +
        scale_x_continuous(expand=c(0, 0)) + scale_y_continuous(expand=c(0, 0))
}

ggmedian <- function(yy, xx) {
    ggplot(data.frame(x=xx, y=yy), aes(x, y)) +
        geom_line() +
        xlab("Temperature") + ylab("Regular response") + theme_bw() +
        scale_x_continuous(expand=c(0, 0)) + scale_y_continuous(expand=c(0, 0))
}

ggcompare <- function(allpoints, yy.base, xx) {
    ggplot(data.frame(x=rep(xx, 2), y=c(allpoints[2,], yy.base), ymin=c(allpoints[1,], yy.base), ymax=c(allpoints[3,], yy.base), group=rep(c('full', 'base'), each=length(xx))), aes(x, y, group=group)) +
        geom_ribbon(aes(ymin=ymin, ymax=ymax), alpha=.5) +
        geom_line(aes(colour=group)) +
        xlab("Temperature") + ylab("Regular response") + theme_bw() +
        scale_colour_discrete(name="Assumptions:", breaks=c('base', 'full'), labels=c("Standard", "Geographic")) + theme(legend.position="bottom")
}

ggmatrix <- function(generate.full, generate.base, covarxs, covarys, covarxlabs, covarylabs, xx) {
    df <- rbind(x=c(), y=c(), ymin=c(), ymax=c(), covarx=c(), covary=c(), group=c())
    for (covarx in covarxs) {
        for (covary in covarys) {
            print(c(covarx, covary))
            allpoints <- generate.full(covarx, covary, xx)
            yy.base <- generate.base(covarx, covary, xx)
            df <- rbind(df, data.frame(x=rep(xx, 2), y=c(allpoints[2,], yy.base), ymin=c(allpoints[1,], yy.base), ymax=c(allpoints[3,], yy.base), group=rep(c('full', 'base'), each=length(xx)), covarx=covarxlabs[covarxs == covarx], covary=covarylabs[covarys == covary]))
        }
    }

    ggplot(df, aes(x, y, group=paste(group, covarx, covary))) +
        facet_grid(covary ~ covarx) +
        geom_ribbon(aes(ymin=ymin, ymax=ymax), alpha=.5) +
        geom_line(aes(colour=group)) +
        xlab("Temperature") + ylab("Regular response") + theme_bw() +
        scale_colour_discrete(name="Assumptions:", breaks=c('base', 'full'), labels=c("Standard", "Geographic")) + theme(legend.position="bottom") +
        scale_x_continuous(expand=c(0, 0)) + scale_y_continuous(expand=c(0, 0))
}
