##
## plot the "SNP Index" for each seg call on each chromosome, stacked.
##
    
strain = "B73"
minsize = 50

opar = par()
par(mfrow=c(10,1),mar=c(0.4,4,0.4,0.4))
    
xmax = max(seg$start[seg$contig==1])

xlim = c(0,xmax)
ylim = c(0,1)

## one plot per chromosome
for (chr in 1:10) {
    
    plot(seg$start[seg$contig==chr], seg$snpindex.nt[seg$contig==chr],
         pch=1, cex=0.5, col="black", ylab=paste("Chr",chr), xlim=xlim, ylim=ylim,  xaxt='n', xaxs='i')
    points(seg$start[seg$contig==chr], seg$snpindex.t[seg$contig==chr],
           pch=1, cex=0.5, col="red", xaxt='n', xaxs='i')
    
}

par(opar)

