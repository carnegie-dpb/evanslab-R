##
## plot log odds ratio of each SNP on each chromosome, stacked plots
##
## chr   pos ref alt twHet twHom twNC ntHet ntHom ntNC log10OR

plot.vcf.logOR.region = function(chr=1, start=0, end=0, ymin=0, ymax=0, minHetSamples=25,
                                 fitCenter1=0, fitHeight1=0, fitWidth1=0,
                                 fitCenter2=0, fitHeight2=0, fitWidth2=0) {
    
    if (start==0) {
        start = min(vcfcounts$pos[vcfcounts$chr==chr])
    }
    if (end==0) {
        end = max(vcfcounts$pos[vcfcounts$chr==chr])
    }
    xlim = c(start,end)

    pts = vcfcounts$chr==chr & vcfcounts$pos>=start & vcfcounts$pos<=end & is.finite(vcfcounts$log10OR) & (vcfcounts$twHet+vcfcounts$ntHet>=minHetSamples)

    if (ymin==0) {
        ymin = min(vcfcounts$log10OR[pts])
    }
    if (ymax==0) {
        ymax = max(vcfcounts$log10OR[pts])
    }
    ylim = c(ymin,ymax)

    plot(vcfcounts$pos[pts], vcfcounts$log10OR[pts],
         xlim=xlim, ylim=ylim,
         xlab=paste("Chr",chr), ylab="log10(O.R.)",
         pch=1, cex=0.5, col="black", 
         main="Twin Study: GBS")
    
    lines(xlim, c(0,0), col="gray")

    ## Gaussian fit
    if (fitCenter1!=0 && fitHeight1!=0 && fitWidth1!=0) {
        fit = fitHeight1*exp(-(vcfcounts$pos[pts]-fitCenter1)^2/fitWidth1^2)
        legend = paste("Fit center =", format(fitCenter1,scientific=FALSE), "width =", format(fitWidth1,scientific=FALSE,width=9))
        if (fitCenter2!=0 && fitHeight2!=0 && fitWidth2!=0) {
            fit = fit + fitHeight2*exp(-(vcfcounts$pos[pts]-fitCenter2)^2/fitWidth2^2)
            legend = c(legend, paste("Fit center =", format(fitCenter2,scientific=FALSE), "width =", format(fitWidth2,scientific=FALSE,width=9)))
        }
        lines(vcfcounts$pos[pts], fit, col="blue")
        legend(x="topright", legend)
    }


    ## SSR markers
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
