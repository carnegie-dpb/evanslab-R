##
## Mo17 calls: plot log Fisher p value of each call on each chromosome, stacked plots
##
## CM009906.1	306176853
## CM009907.1	243312244
## CM009908.1	240957546
## CM009909.1	249748977
## CM009910.1	220382597
## CM009911.1	168720447
## CM009912.1	178574342
## CM009913.1	182612305
## CM009914.1	164939053
## CM009915.1	149041351
    
opar = par(mfrow=c(10,1))
par(mar=c(0.4,4,0.4,0.4))
    
xmax = max(Mo17$start[Mo17$contig=="CM009906.1"])
ymax = 11

xlim = c(0,xmax)
ylim = c(2,ymax)

yline1 = rep(-log10(1e-2), 2)
yline2 = rep(-log10(1e-6), 2)

min = 50

## one plot per chromosome
for (chr in c("CM009906.1", "CM009907.1", "CM009908.1", "CM009909.1", "CM009910.1", "CM009911.1", "CM009912.1", "CM009913.1", "CM009914.1", "CM009915.1")) {
    plot(Mo17$start[Mo17$contig==chr & Mo17$p<1e-2 & Mo17$size>=min], Mo17$mlog10p[Mo17$contig==chr & Mo17$p<1e-2 & Mo17$size>=min],
         pch=1, cex=0.5, col="black", ylab=chr, xlim=xlim, ylim=ylim,  xaxt='n', xaxs='i')
    
    ## highlight significant p values
    points(Mo17$start[Mo17$contig==chr & Mo17$p<1e-6 & Mo17$size>=min], Mo17$mlog10p[Mo17$contig==chr & Mo17$p<1e-6 & Mo17$size>=min], pch=19, cex=0.5, col="darkred")
    lines(xlim,  yline1, col="gray", lty=2)
    lines(xlim,  yline2, col="gray", lty=2)
    
    ## special cases
    if (chr=="CM009906.1") {
        ## title
        text(1, 0.85*ylim[2], cex=1.4, pos=4, "Mo17 calls: Fisher's Exact Test -log10(p) for ALT/REF counts vs Twin/Non-Twin")
    } else if (chr=="CM009907.1") {
        ## filter text
        text(1, 0.85*ylim[2], cex=1.2, pos=4, paste("Calls with at least",min,"total reads"))
    }
}

par(opar)

