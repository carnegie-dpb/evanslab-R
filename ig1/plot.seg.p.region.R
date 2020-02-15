##
## plot -log10(p) of each SNP on each chromosome, stacked plots, for the WGS segregation data
##   contig start REF ALT   a  b   c   d size          p   mlog10p signif
## 1      1  1385   C   T 161 64 162  66  453 0.08212748 1.0855115  false
## a,b are NON-TWIN REF,ALT; c,d are TWIN REF,ALT

plot.seg.p.region = function(chr=1, start=1, end=0, ymin=0.0, ymax=0.0, minReads=10, minAltReads=10) {

    if (end==0) end = max(seg$start[seg$contig==chr])
    xlim = c(start,end)

    if (ymin==0.0) ymin = min(seg$mlog10p[is.finite(seg$mlog10p) & seg$contig==chr])
    if (ymax==0.0) ymax = max(seg$mlog10p[is.finite(seg$mlog10p) & seg$contig==chr])
    ylim = c(ymin,ymax)

    pts = seg$contig==chr & seg$start>=start & seg$start<=end & seg$size>=minReads & (seg$b+seg$d)>=minAltReads
    sigpts = pts & seg$p<1.0e-2

    plot(seg$start[pts], seg$mlog10p[pts],
         xlim=xlim, ylim=ylim,
         xlab=paste("Chr",chr), ylab="-log10(p)", 
         pch=1, cex=0.5, col="black",
         main="Twin Study: BSA"
         )

    ## highlight significant points
    points(seg$start[sigpts], seg$mlog10p[sigpts], pch=19, cex=0.5, col="blue")
    
    ## significance line
    lines(xlim, c(2,2), col="gray")

    ## extra annotation
    if (chr==3) {
        ## show ig1 location
        lines(c(171822155,171822155), c(0,2), col="darkgreen", lwd=3)
        text(171822000, 0, "B73 LBD6 (ig1)", pos=1, col="darkgreen", cex=0.7, offset=0.3)
    }

}
