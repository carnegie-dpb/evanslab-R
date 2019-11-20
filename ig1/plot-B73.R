##
## plot log Fisher p value of each B73 call on each chromosome, stacked plots
##
    
opar = par(mfrow=c(10,1))
par(mar=c(0.4,4,0.4,0.4))
    
xmax = max(B73$start[B73$contig==1])
ymax = 10

xlim = c(0,xmax)
ylim = c(2,ymax)

yline1 = rep(-log10(1e-2), 2)
yline2 = rep(-log10(1e-6), 2)

min = 50

## one plot per chromosome
for (chr in 1:10) {
    
    plot(B73$start[B73$contig==chr & B73$p<1e-2 & B73$size>=min], B73$mlog10p[B73$contig==chr & B73$p<1e-2 & B73$size>=min],
         pch=1, cex=0.5, col="black", ylab=paste("Chr",chr), xlim=xlim, ylim=ylim,  xaxt='n', xaxs='i')
    
    ## highlight highly significant p values
    points(B73$start[B73$contig==chr & B73$p<1e-6 & B73$size>=min], B73$mlog10p[B73$contig==chr & B73$p<1e-6 & B73$size>=min], pch=19, cex=0.5, col="darkred")
    lines(xlim,  yline1, col="gray", lty=2)
    lines(xlim,  yline2, col="gray", lty=2)
    
    ## special cases
    if (chr==1) {
        ## title
        text(1, 0.85*ylim[2], cex=1.5, pos=4, "B73 calls: Fisher's Exact Test -log10(p) for ALT/REF counts vs Twin/Non-Twin")
    } else if (chr==2) {
        ## filter text
        text(1, 0.85*ylim[2], cex=1.2, pos=4, paste("Calls with at least",min,"total reads"))
    } else if (chr==3) {
        ## show B73 location
        lines(c(171820386,171820386), c(0,ymax), col="black")
        lines(c(171823924,171823924), c(0,ymax), col="black")
        text(174500000, ymax, "B73 LBD6 (ig1)", pos=1, col="black")
    } else if (chr==8) {
        ## marker umc1984
        points(c(81869662,81869383), c(ymax-1,ymax-1), col="darkorange", pch=2) 
        text(81869662, ymax-1, "umc1984 (BE050190)", pos=4, col="darkorange")
    } else if (chr==9) {
        ## marker umc1586
        points(c(25250075,25250260), c(ymax-1,ymax-1), col="darkorange", pch=2) 
        text(25250075, ymax-1, "umc1586 (AW053174)", pos=2, col="darkorange")
        ## marker umc2337
        points(c(26675959,26675253), c(ymax-1,ymax-1), col="darkorange", pch=2) 
        text(26675959, ymax-1, "umc2337 (BM339026)", pos=4, col="darkorange")
    }
    
}

par(opar)

