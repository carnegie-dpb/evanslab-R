##
## plot log odds ratio of each SNP on each chromosome, stacked plots
##
## chr   pos ref alt twHet twHom twNC ntHet ntHom ntNC log10OR

plot.vcf.logOR = function(minHetSamples=25) {
    
    opar = par(mfrow=c(10,1))
    par(mar=c(0.4,4,0.4,0.4))
    
    xmax = 305861025
    xmax = max(vcfcounts$pos[vcfcounts$chr==1])
    ymax = 0.5
    
    xlim = c(0,xmax)
    ylim = c(-ymax,ymax)
    
    ## one plot per chromosome
    for (chr in 1:10) {
        pts = vcfcounts$chr==chr & (vcfcounts$twHet+vcfcounts$ntHet)>=minHetSamples
        plot(vcfcounts$pos[pts], vcfcounts$log10OR[pts], pch=1, cex=0.5, col="black", ylab=paste("Chr",chr), xlim=xlim, ylim=ylim, xaxt='n', xaxs='i')
        lines(xlim, c(0,0), col="gray")
    }

    par(opar)
}
