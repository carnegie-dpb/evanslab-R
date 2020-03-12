## plot p-value of each SNP for a given chromosome and region
##
## contig    pos ref alt TWRF TWRR TWAF TWAR NTRF NTRR NTAF NTAR TWRef TWHet  TWHom NTRef NTHet NTHom        p      mlog10p  OR log10OR
##
plot.gbs.p.region = function(chr=1, start=0, end=0, ymin=0, ymax=0, minHetSamples=25,
                             gene=NULL, showGenes=FALSE,
                             fitCenter1=0, fitHeight1=0, fitWidth1=0,
                             fitCenter2=0, fitHeight2=0, fitWidth2=0) {
    if (!is.null(gene)) {
        start = genes[gene,]$start
        end = genes[gene,]$end
    } else if (end==0) {
        end = max(gbscounts$pos[gbscounts$contig==chr])
    }
    xlim = c(start,end)
    
    pts = gbscounts$contig==chr & gbscounts$pos>=start & gbscounts$pos<=end & is.finite(gbscounts$mlog10p) & (gbscounts$TWHet+gbscounts$NTHet)>=minHetSamples

    if (ymax==0) {
        ymax = max(gbscounts$mlog10p[pts])
    }
    ylim = c(ymin,ymax)

    plot(gbscounts$pos[pts], gbscounts$mlog10p[pts],
         xlim=xlim, ylim=ylim,
         xlab=paste("Chr",chr), ylab="-log10(p)",
         pch=1, cex=0.8, col="black", 
         main="Twin Study: GBS")

    ## basic legend
    legend = c(
        paste("min HET samples =", minHetSamples)
    )
    legend(x="topleft", legend, text.col="black")

    
    ## Gaussian fit
    if (fitCenter1!=0 && fitHeight1!=0 && fitWidth1!=0) {
        fit = fitHeight1*exp(-(gbscounts$pos[pts]-fitCenter1)^2/fitWidth1^2)
        legend = paste("Fit center =", format(fitCenter1,scientific=FALSE), "width =", format(fitWidth1,scientific=FALSE,width=9))
        if (fitCenter2!=0 && fitHeight2!=0 && fitWidth2!=0) {
            fit = fit + fitHeight2*exp(-(gbscounts$pos[pts]-fitCenter2)^2/fitWidth2^2)
            legend = c(legend, paste("Fit center =", format(fitCenter2,scientific=FALSE), "width =", format(fitWidth2,scientific=FALSE,width=9)))
        }
        lines(gbscounts$pos[pts], fit, col="blue")
        legend(x="topright", legend)
    }

    ## extra annotation on chr 3
    if (chr==3) {
        ## show ig1 location
        ig1start = 171820603;
        ig1end = 171823900;
        lines(rep((ig1start+ig1end)/2,2), c(ymin,ymax), col="darkgreen", lwd=3)
        text((ig1start+ig1end)/2, ymin, "ig1", pos=1, offset=0.4, cex=0.7, col="darkgreen")
        ## show genetic markers
        source("load-markers.R")
        markers = markers[markers$chr==chr,]
        for (i in 1:nrow(markers)) {
            lines(c(markers$start[i],markers$start[i]), c(ymin,ymax), col="blue", lty=3, lwd=1)
            text(markers$start[i], ymin, markers$marker[i], col="blue", pos=3, cex=0.5)
            text(markers$start[i], ymin-0.05*(ymax-ymin), paste(markers$coord[i],"cM"), col="blue", pos=3, cex=0.5)
        }
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
