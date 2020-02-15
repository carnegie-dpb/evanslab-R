library("zoo")

## plot log odds ratio of each SNP on each chromosome, stacked plots, for the WGS segregation data
##   contig start REF ALT   a  b   c   d size          p   mlog10p signif
## 1      1  1385   C   T 161 64 162  66  453 0.08212748 1.0855115  false
## a,b are NON-TWIN REF,ALT; c,d are TWIN REF,ALT

plot.both.logOR.region = function(chr=1, start=1, end=0, ymin=0, ymax=0, minAltReads=20, minHetSamples=25, nAvg=0) {

    if (end==0) {
        end = max(vcfcounts$pos[vcfcounts$chr==chr])
    }
    xlim = c(start,end)

    vcf.pts = vcfcounts$chr==chr & is.finite(vcfcounts$log10OR) & vcfcounts$pos>=start & vcfcounts$pos<=end & (vcfcounts$twHet+vcfcounts$ntHet>=minHetSamples)
    seg.pts = seg$contig==chr & is.finite(seg$log10OR) & seg$start>=start & seg$start<=end & (seg$b+seg$d)>=minAltReads


    if (ymin==0) ymin = min(c(vcfcounts$log10OR[vcf.pts], seg$log10OR[seg.pts]))
    if (ymax==0) ymax = max(c(vcfcounts$log10OR[vcf.pts], seg$log10OR[seg.pts]))
    ylim = c(ymin,ymax)

    ## index seg$start, etc. contiguously from 1 to enable averaging below
    starts = seg$start[seg.pts]
    a = seg$a[seg.pts]
    b = seg$b[seg.pts]
    c = seg$c[seg.pts]
    d = seg$d[seg.pts]
    x = c()
    y = c()

    if (nAvg>1) {
        for (i in seq(1,length(starts),nAvg)) {
            range = seq(i,i+nAvg-1)
            x = c(x, mean(starts[range]))
            a.sum = sum(a[range])
            b.sum = sum(b[range])
            c.sum = sum(c[range])
            d.sum = sum(d[range])
            y = c(y, log10((d.sum/c.sum)/(b.sum/a.sum)))
        }
    } else {
        x = starts
        y = seg$log10OR[seg.pts]
    }

    plot(x, y, pch=19, cex=0.3, col="blue",
           xlab=paste("Chr",chr), ylab="log10(O.R.)", xlim=xlim, ylim=ylim, main="GBS (red) and WGS (blue)")
    
    points(vcfcounts$pos[vcf.pts], vcfcounts$log10OR[vcf.pts], pch=19, cex=0.5, col="red")

    lines(xlim, c(0,0), col="gray")
}

