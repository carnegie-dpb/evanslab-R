##
## plot log odds ratio of each SNP on each chromosome, stacked plots, for the WGS segregation data
##   contig start REF ALT   a  b   c   d size          p   mlog10p signif
## 1      1  1385   C   T 161 64 162  66  453 0.08212748 1.0855115  false
## a,b are NON-TWIN REF,ALT; c,d are TWIN REF,ALT

plot.seg.logOR.region = function(chr=1, start=1, end=0, ymin=0.0, ymax=0.0, minReads=10, minAltReads=10,
                                 fitCenter1=0, fitHeight1=0, fitWidth1=0,
                                 fitCenter2=0, fitHeight2=0, fitWidth2=0) {

    if (end==0) end = max(seg$start[seg$contig==chr])
    xlim = c(start,end)

    if (ymin==0.0) ymin = min(seg$log10OR[is.finite(seg$log10OR) & seg$contig==chr])
    if (ymax==0.0) ymax = max(seg$log10OR[is.finite(seg$log10OR) & seg$contig==chr])
    ylim = c(ymin,ymax)

    pts = seg$contig==chr & seg$start>=start & seg$start<=end & seg$size>=minReads & (seg$b+seg$d)>=minAltReads
    sigpts = pts & seg$p<1.0e-2

    plot(seg$start[pts], seg$log10OR[pts],
         xlim=xlim, ylim=ylim,
         xlab=paste("Chr",chr), ylab="log10(O.R.)", 
         pch=1, cex=0.5, col="black",
         main="Twin Study: BSA"
         )

    ## highlight significant points
    points(seg$start[sigpts], seg$log10OR[sigpts], pch=19, cex=0.5, col="blue")
    
    ## zero line
    lines(xlim, c(0,0), col="gray")
    
    ## Gaussian fit(s)
    if (fitCenter1!=0 && fitHeight1!=0 && fitWidth1!=0) {
        fit = fitHeight1*exp(-(seg$start[pts]-fitCenter1)^2/fitWidth1^2)
        legend = paste("Fit center =", chr, ":", format(fitCenter1,scientific=FALSE), "width =", format(fitWidth1,scientific=FALSE,width=9))
        if (fitCenter2!=0 && fitHeight2!=0 && fitWidth2!=0) {
            fit = fit + fitHeight2*exp(-(seg$start[pts]-fitCenter2)^2/fitWidth2^2)
            legend = c(legend, paste("Fit center =", chr, ":", format(fitCenter2,scientific=FALSE), "width =", format(fitWidth2,scientific=FALSE,width=9)))
        }
        lines(seg$start[pts], fit, col="darkgreen", lwd=3)
        legend(x="topleft", legend, text.col="darkgreen")
    }
    
    ## extra annotation
    if (chr==3) {
        ## show ig1 location
        lines(c(171822155,171822155), c(-0.5,0), col="darkgreen", lwd=3)
        text(171822000, -0.5, "B73 LBD6 (ig1)", pos=1, col="darkgreen", cex=0.7, offset=0.3)
    }
}
