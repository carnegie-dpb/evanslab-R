##
## plot the "SNP Index" for each seg call on each chromosome, stacked.
##
    
opar = par()
par(mfrow=c(10,1),mar=c(0.4,4,0.4,0.4))

strain = "B73"
minsize = 20
    
xmax = max(seg$start[seg$contig==1])

xlim = c(0,xmax)
ylim = c(-1,1)

plotrange = (seg$a+seg$b)>minsize & (seg$c+seg$d)>minsize
plotseg = seg[plotrange,]

## one plot per chromosome
for (chr in 1:10) {
    plot(plotseg$start[plotseg$contig==chr],    plotseg$snpindex.t[plotseg$contig==chr],
         pch=1, cex=0.2, col="darkred", ylab=paste("Chr",chr), xlim=xlim, ylim=ylim,  xaxt='n', xaxs='i')
    points(plotseg$start[plotseg$contig==chr], -plotseg$snpindex.nt[plotseg$contig==chr],
           pch=1, cex=0.2, col="darkblue")
}

par(opar)

