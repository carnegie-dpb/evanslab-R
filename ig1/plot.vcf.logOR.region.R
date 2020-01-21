##
## plot log odds ratio of each SNP on each chromosome, stacked plots
##
## chr   pos ref alt twHet twHom twNC ntHet ntHom ntNC log10OR

plot.vcf.logOR.region = function(chr=1, start=1, end=1000000000, ymax=0.2, minHetSamples=25) {
    
    xlim = c(start,end)
    ylim = c(-ymax,ymax)
    
    ## one plot per chromosome
    pts = vcfcounts$chr==chr & vcfcounts$pos>=start & vcfcounts$pos<=end & (vcfcounts$twHet+vcfcounts$ntHet>=minHetSamples)

    plot(vcfcounts$pos[pts], vcfcounts$log10OR[pts], pch=1, cex=0.5, col="black", xlab=paste("Chr",chr), ylab="-log10(O.R.)", main="GBS", xlim=xlim, ylim=ylim)
    lines(xlim, c(0,0), col="gray")
}
