##
## plot -log10(p) of each SNP on each chromosome, stacked plots, for the WGS segregation data
##   contig start REF ALT   a  b   c   d size          p   mlog10p signif
## 1      1  1385   C   T 161 64 162  66  453 0.08212748 1.0855115  false
## a,b are NON-TWIN REF,ALT; c,d are TWIN REF,ALT

plot.seg.p.region = function(chr=1, start=1, end=0, ymin=0, ymax=0, minReads=20, minAltReads=10) {

    if (end==0) end = max(seg$start[seg$contig==chr])
    xlim = c(start,end)

    pts = seg$contig==chr & seg$start>=start & seg$start<=end & seg$size>=minReads & (seg$b+seg$d)>=minAltReads & is.finite(seg$mlog10p)

    if (ymin==0) {
        ymin = min(seg$mlog10p[pts])
    }
    if (ymax==0.0) {
        ymax = max(seg$mlog10p[pts])
    }
    ylim = c(ymin,ymax)

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

    ## basic legend
    legend = c(paste("min. reads =", minReads),
               paste("min. ALT reads =", minAltReads))
    legend(x="topleft", legend, text.col="black")

    ## extra annotation
    if (chr==3) {
        ## show ig1 location
        lines(c(171822155,171822155), c(ymin,ymax*0.9), col="darkgreen", lwd=3)
        text(171822000, ymax, "B73 LBD6 (ig1)", pos=1, col="darkgreen", cex=1.0)
        ## Name		gene:Zm00001d043330
        ## Description	maternal effect embryo arrest 22
        ## Position	3:196835213..196846471 (- strand)
        lines(c(196835213,196846471), c(ymin,ymin), col="blue", lty=1, lwd=2)
        text((196835213+196846471)/2, ymin, "Zm00001d043330", pos=1, offset=0.4, cex=0.7, col="blue")
        text(196835213, ymin, "<", pos=4, offset=0, col="blue", cex=1)
        text(196846471, ymin, "<", pos=2, offset=0, col="blue", cex=1)
        ## Name         gene:Zm00001d043333
        ## Position     3:196851706..196853791 (+ strand)
        ## Length       2,086 bp
        lines(c(196851706,196853791), c(ymin,ymin), col="blue", lty=1, lwd=2)
        text((196851706+196853791)/2, ymin, "Zm00001d043333", pos=1, offset=0.4, cex=0.7, col="blue")
        text(196851706, ymin, ">", pos=4, offset=0, col="blue", cex=1)
        text(196853791, ymin, ">", pos=2, offset=0, col="blue", cex=1)
    }
}
