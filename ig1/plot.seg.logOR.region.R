##
## plot log odds ratio of each SNP on each chromosome, stacked plots, for the WGS segregation data
##   contig start REF ALT   a  b   c   d size          p   mlog10p signif
## 1      1  1385   C   T 161 64 162  66  453 0.08212748 1.0855115  false
## a,b are NON-TWIN REF,ALT; c,d are TWIN REF,ALT

plot.seg.logOR.region = function(chr=1, start=1, end=1e8, ymax=0.2, minAltReads=10) {
    
    xlim = c(start,end)
    ylim = c(-ymax,ymax)

    pts = seg$contig==chr & seg$start>=start & seg$start<=end & seg$b>=minAltReads & seg$d>=minAltReads

    plot(seg$start[pts], seg$log10OR[pts], pch=1, cex=0.5, col="black", xlab=paste("Chr",chr), ylab="-log10(O.R.)", xlim=xlim, ylim=ylim)
    lines(xlim, c(0,0), col="gray")

}
