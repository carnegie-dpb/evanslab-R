##
## plot -log10(p) of each SNP on each chromosome, stacked plots, for the BSA
##
plot.seg.p.region = function(chr=1, start=1, end=0, ymin=0, ymax=0, minAltReads=20, maxLogRatio2=4.0, pSig=5.0e-8, showGenes=FALSE, gene=NULL) {

    if (!is.null(gene)) {
        start = genes[gene,]$start
        end = genes[gene,]$end
    } else if (end==0) {
        end = max(seg$start[seg$contig==chr])
    }
    xlim = c(start,end)

    pts = seg$contig==chr & seg$start>=start & seg$start<=end & is.finite(seg$mlog10p) & (seg$alt1+seg$alt2)>=minAltReads

    ## require legit HET calls on twins
    pts = pts & abs(seg$logRatio2)<maxLogRatio2

    if (ymax==0) {
        ymax = max(seg$mlog10p[pts])
    }
    ylim = c(ymin,ymax)

    sigpts = pts & seg$p<pSig

    plot(seg$start[pts], seg$mlog10p[pts],
         xlim=xlim, ylim=ylim,
         xlab=paste("Chr",chr), ylab="-log10(p)", 
         pch=1, cex=0.5, col="black",
         main="Twin Study: BSA"
         )

    ## highlight significant points
    points(seg$start[sigpts], seg$mlog10p[sigpts], pch=19, cex=0.5, col="blue")
    
    ## basic legend
    legend = c(
        paste("min ALT reads =", minAltReads),
        paste("max |log2(ALT2/REF2)| =", maxLogRatio2)
    )
    legend(x="topleft", legend, text.col="black")

    ## extra annotation
    if (chr==3) {
        ## show ig1 location
        lines(c(171822155,171822155), c(ymin,ymax*0.9), col="darkgreen", lwd=3)
        text(171822000, ymin, "ig1", pos=1, offset=0.4, cex=0.7, col="darkgreen")
    }

    if (!is.null(gene) || showGenes) {
        ##                seqid  source type  start    end
        ## Zm00001d027230     1 gramene gene  44289  49837
        inGenes = genes[genes$seqid==chr & genes$end>=start & genes$start<=end,]
        for (i in 1:nrow(inGenes)) {
            lines(c(inGenes$start[i],inGenes$end[i]), c(ymin,ymin), col="blue", lwd=4)
            text((inGenes$start[i]+inGenes$end[i])/2, ymin, rownames(inGenes)[i], pos=1, offset=0.4, cex=0.8, col="blue")
            if (inGenes$strand[i]=="+") {
                text(inGenes$end[i], ymin, ">", pos=4, offset=-0.3, cex=1.5, col="blue")
            } else {
                text(inGenes$start[i], ymin, "<", pos=2, offset=-0.3, cex=1.5, col="blue")
            }
            geneCDS = cds[cds$gene==rownames(inGenes)[i],]
            for (j in 1:nrow(geneCDS)) {
                lines(c(geneCDS$start[j],geneCDS$end[j]), c(ymin,ymin), col="orange", lwd=4)
            }
        }
    }
}
