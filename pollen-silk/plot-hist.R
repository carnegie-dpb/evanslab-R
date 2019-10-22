## various ways to plot a histogram of counts
## requires mergeVCF.R beforehand!

library("MASS")

## overplot a log-normal distribution REF/REF fit on a hist plot for pollen ratio estimate
hist(log(merged.bigref$REFtot.m/merged.bigref$REFtot.r),
     breaks=30, freq=FALSE,
     xlab=paste("log(",mixSampleName,"REF /",refSampleName,"REF)"),
     main="B73 map REF reads ratio (bigref loci)")
med = median(merged.bigref$REFtot.m/merged.bigref$REFtot.r)
fit = fitdistr(merged.bigref$REFtot.m/merged.bigref$REFtot.r, "lognormal")


## ## overplot a log-normal distribution TOT/TOT fit on a hist plot for total coverage ratio estimate
## hist(log(merged.strong$TOT.m/merged.strong$TOT.r),
##      breaks=30, freq=FALSE,
##      ylim=c(0,.2),
##      xlab=paste("log(",mixSampleName,"TOT /",refSampleName,"TOT)"),
##      main="B73 map TOT reads ratio")
## med = median(merged.strong$TOT.m/merged.strong$TOT.r)
## fit = fitdistr(merged.strong$TOT.m/merged.strong$TOT.r, "lognormal")


## ## overplot a log-normal distribution REF.m / ALT.m and fit on a hist plot to estimate pollen dilution = ratio / (1 + ratio)
## hist(log(merged.bigref$REFtot.m/merged.bigref$ALTtot.m),
##      breaks=30, freq=FALSE,
##      xlab=paste("log(",mixSampleName,"REF /",mixSampleName," ALT)"),
##      main="B73 map MIX REF / MIX ALT reads")
## med = median(merged.bigref$REFtot.m/merged.bigref$ALTtot.m)
## fit = fitdistr(merged.bigref$REFtot.m/merged.bigref$ALTtot.m, "lognormal")

## annotation
lines(log(c(med,med)), c(0,1), col="darkgreen", lwd=2)
meanlog = fit[[1]][1]
sdlog = fit[[1]][2]
meanlog.err = fit[[2]][1]
sdlog.err = fit[[2]][2]
mean = exp(meanlog)
mean.err = mean - exp(meanlog-meanlog.err)

x = c(-100:100/10)
lines(x, exp(-((x-meanlog)/sdlog)^2/2)/(sdlog*sqrt(2*pi)))
text(log(med), 1.0/(sdlog*sqrt(2*pi)), paste("median ratio =", round(med,3)), col="darkgreen", pos=2, offset=2)
text(log(med), 1.0/(sdlog*sqrt(2*pi)), paste("fit mean =", round(mean,3), "+/-", round(mean.err,3)), col="darkgreen", pos=4, offset=2)
