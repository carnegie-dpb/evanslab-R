library("zoo")

## plot log p each SNP on each chromosome, stacked plots, for the WGS segregation data
##   contig start REF ALT   a  b   c   d size          p   mlog10p signif
## 1      1  1385   C   T 161 64 162  66  453 0.08212748 1.0855115  false
## a,b are NON-TWIN REF,ALT; c,d are TWIN REF,ALT

plot.both.p.region = function(chr=1, start=1, end=0, ymin=0, ymax=0, minAltReads=20, minHetSamples=25) {

    if (end==0) {
        end = max(vcfcounts$pos[vcfcounts$chr==chr])
    }
    xlim = c(start,end)

    vcf.pts = vcfcounts$chr==chr & is.finite(vcfcounts$mlog10p) & vcfcounts$pos>=start & vcfcounts$pos<=end & (vcfcounts$twHet+vcfcounts$ntHet>=minHetSamples)
    seg.pts = seg$contig==chr & is.finite(seg$mlog10p) & seg$start>=start & seg$start<=end & (seg$b+seg$d)>=minAltReads

    if (ymin==0) ymin = min(c(vcfcounts$mlog10p[vcf.pts], seg$mlog10p[seg.pts]))
    if (ymax==0) ymax = max(c(vcfcounts$mlog10p[vcf.pts], seg$mlog10p[seg.pts]))
    ylim = c(ymin,ymax)

    plot(seg$start[seg.pts], seg$mlog10p[seg.pts], pch=19, cex=0.3, col="blue",
         xlab=paste("Chr",chr), ylab="-log10(p)", xlim=xlim, ylim=ylim, main="GBS (red) and WGS (blue)")
    
    points(vcfcounts$pos[vcf.pts], vcfcounts$mlog10p[vcf.pts], pch=19, cex=0.5, col="red")
    
    lines(xlim, c(0,0), col="gray")
}

