##
## plot log odds ratio of each SNP on each chromosome, stacked plots, for the WGS segregation data
##   contig start REF ALT   a  b   c   d size          p   mlog10p signif
## 1      1  1385   C   T 161 64 162  66  453 0.08212748 1.0855115  false
##
plot.seg.or = function(minReads=4, minHetSamples=45) {
    
    opar = par(mfrow=c(10,1))
    par(mar=c(0.4,4,0.4,0.4))
    
    xmax = 305861025
    xmax = max(seg$start[seg$contig==1])
    ymax = 2
    
    xlim = c(0,xmax)
    ylim = c(-ymax,ymax)
    
    yline1 = rep(0, 2)
    yline2 = rep(log10(3), 2)
    
    ## sample counts
    nsamples = max(seg$Ahet+seg$Ahom+seg$Bhet+seg$Bhom)
    nAsamples = max(seg$Ahom+seg$Ahet) # hopefully the number of non-twin samples
    nBsamples = max(seg$Bhom+seg$Bhet) # hopefully the number of non-twin samples

    ## require a minimum number of Het samples
    filter = (seg$Ahet+seg$Bhet)>=minHetSamples

    ## require a minimum average number of reads per CALLED sample
    if (minReads>0) {
        filter = filter & ((seg$Aref+seg$Aalt+seg$Bref+seg$Balt)/(seg$Ahom+seg$Ahet+seg$Bhom+seg$Bhet)>=minReads)
    }

    ## Odds Ratio
    AhetOdds = seg$Ahet/(nAsamples-seg$Ahet)
    BhetOdds = seg$Bhet/(nBsamples-seg$Bhet)
    OR = BhetOdds/AhetOdds

    ## one plot per chromosome
    for (chr in 1:10) {

        plot(seg$start[seg$contig==chr & filter], log10(OR[seg$contig==chr & filter]), pch=1, col="black", ylab=paste("Chr",chr), xlim=xlim, ylim=ylim,  xaxt='n', xaxs='i')
        
        ## highlight significant p values
        points(seg$start[seg$contig==chr & filter & OR>3], log10(OR[seg$contig==chr & filter & OR>3]), pch=19, col="darkred")
        lines(xlim,  yline1, col="gray", lty=2)
        lines(xlim,  yline2, col="gray", lty=2)
        
        ## special cases
        if (chr==1) {
            ## title
            text(1, 0.80*ylim[2], cex=1.5, pos=4, "log10(Odds Ratio) Twin Het:Twin NC/Non-Twin Het:Non-Twin NC")
        } else if (chr==2) {
            ## filter text
            subtitle = paste("Plotting filter: require at least",minHetSamples,"Het samples out of",nsamples)
            if (minReads>0) subtitle = paste(subtitle,"; loci must have at least",minReads,"average reads per called sample.")
            text(1, 0.80*ylim[2], cex=1.5, pos=4, subtitle)
        } else if (chr==3) {
            ## show ig1 location
            lines(c(171820386,171820386), c(0,ymax), col="black")
            lines(c(171823924,171823924), c(0,ymax), col="black")
            text(174500000, ymax, "B73 LBD6 (ig1)", pos=1, col="black")
        } else if (chr==8) {
            ## marker umc1984
            points(c(81869662,81869383), c(ymax-1,ymax-1), col="darkorange", pch=2) 
            text(81869662, 0.80*ymax, "umc1984 (BE050190)", pos=4, col="darkorange")
        } else if (chr==9) {
            ## marker umc1586
            points(c(25250075,25250260), c(0.80*ymax,0.80*ymax), col="darkorange", pch=2) 
            text(25250075, 0.80*ymax, "umc1586 (AW053174)", pos=2, col="darkorange")
            ## marker umc2337
            points(c(26675959,26675253), c(0.80*ymax,0.80*ymax), col="darkorange", pch=2) 
            text(26675959, 0.80*ymax, "umc2337 (BM339026)", pos=4, col="darkorange")
        }
        
    }

    par(opar)

}
