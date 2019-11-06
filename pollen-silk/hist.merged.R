## various ways to plot a histogram of counts
## run load-merged.R beforehand!

library("MASS")

hist.merged = function(altAllele=FALSE) {

    if (altAllele) {
        allele = "ALT"
        numerator =   merged$ALTtot.m[bigalt]
        denominator = merged$ALTtot.a[bigalt]
        xLabel = paste("log(",mixSampleName,"ALT /",altSampleName,"ALT)")
        mainTitle = paste(reference,"map ALT reads ratio (bigalt loci)")
    } else {
        allele = "REF"
        numerator =   merged$REFtot.m[bigref]
        denominator = merged$REFtot.r[bigref]
        xLabel = paste("log(",mixSampleName,"REF /",refSampleName,"REF)")
        mainTitle = paste(reference,"map REF reads ratio (bigref loci)")
    }

    ## overplot a log-normal distribution REF/REF fit on a hist plot for REF ratio estimate
    hist(log(numerator/denominator),
         breaks=31, freq=FALSE,
         xlab=xLabel,
         main=mainTitle
         )

    ## fits
    med = median(numerator/denominator)
    fit = fitdistr(numerator/denominator, "lognormal")

    ## annotation
    meanlog = as.numeric(fit[[1]][1])
    sdlog = as.numeric(fit[[1]][2])
    meanlog.err = as.numeric(fit[[2]][1])
    sdlog.err = as.numeric(fit[[2]][2])
    mean = exp(meanlog)
    mean.err = mean - exp(meanlog-meanlog.err)

    lines(log(c(med,med)), c(0,1), col="darkblue", lwd=2)
    lines(c(meanlog,meanlog), c(0,1), col="darkgreen", lwd=2)

    x = c(-100:100/10)
    lines(x, exp(-((x-meanlog)/sdlog)^2/2)/(sdlog*sqrt(2*pi)))
    text(log(med), par()$yaxp[2], paste("median ratio =", round(med,3)), col="darkblue", pos=2, offset=1)
    text(meanlog, par()$yaxp[2], paste("fit mean =", round(mean,5), "+/-", round(mean.err,5)), col="darkgreen", pos=4, offset=1)
}
