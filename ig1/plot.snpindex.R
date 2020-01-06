##
## plot the "SNP Index" over a given region
##

plot.snpindex = function(chr, start, end, geneIds=FALSE, boxcar=5, diff=FALSE) {
    
    strain = "B73"
    minsize = 100
    xlim = c(start,end)

    plotrange = seg$contig==chr & seg$start>=start & seg$start<=end
    plotseg = seg[plotrange,]

    box.nt = rep(0.5,boxcar)
    box.t = rep(0.5,boxcar)
    imin = boxcar + 1
    imax = nrow(plotseg) - boxcar
    ## add up the individual counts and divide those rather than averaging snpindex values
    for (i in imin:imax) {
        ilims = (i-boxcar):(i+boxcar)
        a = sum(plotseg$a[ilims])
        b = sum(plotseg$b[ilims])
        c = sum(plotseg$c[ilims])
        d = sum(plotseg$d[ilims])
        if ((a+b)>minsize && (c+d)>minsize) {
            box.t = c(box.t, b/(a+b))
            box.nt = c(box.nt, d/(c+d))
        } else {
            box.nt = c(box.nt, 0.5)
            box.t = c(box.t, 0.5)
        }
    }
    box.nt = c(box.nt, rep(0.5,boxcar))
    box.t = c(box.t, rep(0.5,boxcar))

    if (diff) {
        ylim = c(-1,1)
        plot(plotseg$start, plotseg$snpindex.t-plotseg$snpindex.nt,
             pch=1, cex=0.2, col="black", xlab="Position", ylab=paste("Chr",chr), xlim=xlim, ylim=ylim)
        lines(plotseg$start, box.t-box.nt, col="blue")
        title(main=paste(strain,"calls: SNP-index Twin MINUS Non-twin\nCalls with at least",minsize,"total reads"), cex.main=0.9)
    } else {
        ylim = c(0,1)
        plot(plotseg$start, plotseg$snpindex.nt, pch=1, cex=0.2, col="blue",
             xlab="Position", ylab=paste("Chr",chr), xlim=xlim, ylim=ylim)
        points(plotseg$start, plotseg$snpindex.t, pch=1, cex=0.2, col="red")
        lines(plotseg$start, box.nt, col="darkblue", lwd=2)
        lines(plotseg$start, box.t, col="darkred", lwd=2)
        title(main=paste(strain,"calls: SNP-index for Twin (red) and Non-twin (blue)\nCalls with at least",minsize,"total reads"), cex.main=0.9)
    }
    
}
