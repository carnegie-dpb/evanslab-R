## plot log odds ratio of each SNP for a given chromosome and region
##
## contig   pos ref alt TWRF TWRR TWAF TWAR NTRF NTRR NTAF NTAR TWRef TWHet TWHom NTRef NTHet NTHom        p   mlog10p  OR log10OR

plot.gbs.p.region = function(chr=1, start=0, end=0, ymin=0, ymax=0, minHetSamples=25,
                              fitCenter1=0, fitHeight1=0, fitWidth1=0,
                              fitCenter2=0, fitHeight2=0, fitWidth2=0) {
    
    if (start==0) {
        start = min(gbscounts$pos[gbscounts$contig==chr])
    }
    if (end==0) {
        end = max(gbscounts$pos[gbscounts$contig==chr])
    }
    xlim = c(start,end)

    pts = gbscounts$contig==chr & gbscounts$pos>=start & gbscounts$pos<=end & is.finite(gbscounts$mlog10p) & (gbscounts$TWHet+gbscounts$NTHet)>=minHetSamples

    if (ymin==0) {
        ymin = min(gbscounts$mlog10p[pts])
    }
    if (ymax==0) {
        ymax = max(gbscounts$mlog10p[pts])
    }
    ylim = c(ymin,ymax)

    plot(gbscounts$pos[pts], gbscounts$mlog10p[pts],
         xlim=xlim, ylim=ylim,
         xlab=paste("Chr",chr), ylab="-log10(p)",
         pch=1, cex=1.0, col="black", 
         main="Twin Study: GBS")
    
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
        ## text(start, ymax, "Markers: SSR B73 x Mo17 Davis 3 map", pos=4)
        ## ig1 marker
        ## 3 gramene gene 171820603 171823900 ID=gene:Zm00001d042560;Name=LBD6,ig1
        lines(c(171820603,171820603), c(ymin,ymax), col="red", lty=2, lwd=1)
        text(171823900, ymax, "ig1", col="red", pos=4)
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
