##
## plot the "SNP Index" over a given region
##

plot.snpindex = function(chr, start, end, geneIds=FALSE, diff=FALSE, twinOnly=FALSE, nontwinOnly=FALSE, average=0) {
    
    strain = "B73"
    minsize = 20
    xlim = c(start,end)

    plotrange = seg$contig==chr & seg$start>=start & seg$start<=end & (seg$a+seg$b)>minsize & (seg$c+seg$d)>minsize
    plotseg = seg[plotrange,]

    ## add up ALL the individual counts and divide those rather than averaging snpindex values
    if (average>0) {
        fullrange = seg$contig==chr & seg$start>=start & seg$start<=end
        fullseg = seg[fullrange,]
        ave.nt = c()
        ave.t = c()
        avestart = c()
        for (i in seq(1, nrow(fullseg), by=average)) {
            irange = i:i+average-1
            if (diff || nontwinOnly || (!twinOnly && !nontwinOnly)) {
                a = sum(fullseg$a[irange])
                b = sum(fullseg$b[irange])
                ave.nt = c(ave.nt, b/(a+b), b/(a+b))
            }
            if (diff || twinOnly || (!twinOnly && !nontwinOnly)) {
                c = sum(fullseg$c[irange])
                d = sum(fullseg$d[irange])
                ave.t = c(ave.t, d/(c+d), d/(c+d))
            }
            avestart = c(avestart, fullseg$start[i], fullseg$start[i+average-1])
        }
    }

    if (diff) {
        ylim = c(-1,1)
        plot(plotseg$start, plotseg$snpindex.t-plotseg$snpindex.nt,
             pch=1, cex=0.2, col="green", ylab="SNP-index Difference", xlab=paste("Chr",chr,"position"), xlim=xlim, ylim=ylim,
             main=paste(strain,"calls: SNP-index Twin MINUS Non-twin\nCalls with at least",minsize,"total reads"), cex.main=0.9)
        if (average>0) {
            lines(avestart, ave.t-ave.nt, col="darkgreen")
        }
    } else if (twinOnly) {
        ylim = c(0,1)
        plot(plotseg$start, plotseg$snpindex.t, pch=1, cex=0.2, col="red",
             ylab="SNP-index", xlab=paste("Chr",chr,"position"), xlim=xlim, ylim=ylim,
             main=paste(strain,"calls: SNP-index for Twin calls with at least",minsize,"reads"), cex.main=0.9)
        if (average>0) {
            lines(avestart, ave.t, col="darkred", lwd=2)
        }
    } else if (nontwinOnly) {
        ylim = c(0,1)
        plot(plotseg$start, plotseg$snpindex.nt, pch=1, cex=0.2, col="blue",
             ylab="SNP-index", xlab=paste("Chr",chr,"position"), xlim=xlim, ylim=ylim,
             main=paste(strain,"calls: SNP-index for Non-twin calls with at least",minsize,"reads"), cex.main=0.9)
        if (average>0) {
            lines(avestart, ave.nt, col="darkblue", lwd=2)
        }
    } else {
        ylim = c(0,1)
        plot(plotseg$start, plotseg$snpindex.t, pch=1, cex=0.2, col="red",
             ylab="SNP-index", xlab=paste("Chr",chr,"position"), xlim=xlim, ylim=ylim,
             main=paste(strain,"calls: SNP-index for Twin (red) and Non-twin (blue)\nCalls with at least",minsize,"reads per sample type"), cex.main=0.9)
        points(plotseg$start, plotseg$snpindex.nt, pch=1, cex=0.2, col="blue")
        if (average>0) {
            lines(avestart, ave.t, col="darkred", lwd=2)
            lines(avestart, ave.nt, col="darkblue", lwd=2)
        }
    }

    genes = gff[gff$type=="gene" & gff$seqid==chr & gff$start>=start & gff$start<=end,]
    if (nrow(genes)>0) {
        for (i in 1:nrow(genes)) {
            lines(c(genes$start[i],genes$end[i]), c(ylim[1],ylim[1]), col="blue", lwd=2)
            if (geneIds) {
                geneId = sub("ID=gene:", "", strsplit(genes$attributes[i],";")[[1]][1])
                text(genes$start[i], ylim[1]+0.05, geneId, cex=0.8, pos=4, col="blue", offset=0)
            }
        }
    }
}
