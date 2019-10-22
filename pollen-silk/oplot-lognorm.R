## overplot a log-normal distribution fit (presumably to a hist plot)

fit = fitdistr(merged.strong$REFtot.m/merged.strong$REFtot.r, "lognormal")
meanlog = fit[[1]][1]
sdlog = fit[[1]][2]

x = c(-100:100/10)
lines(x, exp(-((x-meanlog)/sdlog)^2/2)/(sdlog*sqrt(2*pi)))
