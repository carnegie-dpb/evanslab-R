## plot log odds ratio of each SNP for a given chromosome and region
##
## contig   pos ref alt TWRF TWRR TWAF TWAR NTRF NTRR NTAF NTAR TWRef TWHet TWHom NTRef NTHet NTHom        p   mlog10p  OR log10OR

plot.gbs.or.region = function(chr=1, start=0, end=0, ymin=0, ymax=0, minHetSamples=25,
                              fitCenter1=0, fitHeight1=0, fitWidth1=0,
                              fitCenter2=0, fitHeight2=0, fitWidth2=0) {
    
    if (start==0) {
        start = min(gbscounts$pos[gbscounts$contig==chr])
    }
    if (end==0) {
        end = max(gbscounts$pos[gbscounts$contig==chr])
    }
    xlim = c(start,end)

    pts = gbscounts$contig==chr & gbscounts$pos>=start & gbscounts$pos<=end & is.finite(gbscounts$log10OR) & (gbscounts$TWHet+gbscounts$NTHet>=minHetSamples)

    if (ymin==0) {
        ymin = min(gbscounts$log10OR[pts])
    }
    if (ymax==0) {
        ymax = max(gbscounts$log10OR[pts])
    }
    ylim = c(ymin,ymax)

    plot(gbscounts$pos[pts], gbscounts$log10OR[pts],
         xlim=xlim, ylim=ylim,
         xlab=paste("Chr",chr), ylab="log10(O.R.)",
         pch=1, cex=0.5, col="black", 
         main="Twin Study: GBS")
    
    lines(xlim, c(0,0), col="gray")

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

    ## extra annotation
    if (chr==3) {
        source("load-markers.R")
        markers = markers[markers$chr==chr,]
        for (i in 1:nrow(markers)) {
            lines(c(markers$start[i],markers$start[i]), c(ymin,ymax), col="blue", lty=3, lwd=1)
            text(markers$start[i], ymin, markers$marker[i], col="blue", pos=3, cex=0.5)
            text(markers$start[i], ymin-0.05*(ymax-ymin), paste(markers$coord[i],"cM"), col="blue", pos=3, cex=0.5)
        }
        text(start, ymax, "Markers: SSR B73 x Mo17 Davis 3 map", pos=4)
        ## ig1 marker
        ## 3 gramene gene 171820603 171823900 ID=gene:Zm00001d042560;Name=LBD6,ig1
        lines(c(171820603,171820603), c(ymin,ymax), col="red", lty=2, lwd=1)
        text(171823900, ymax, "ig1", col="red", pos=4)
    }
}
