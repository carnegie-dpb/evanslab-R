##
## plot log Fisher p value of each SNP on each chromosome, stacked plots
##
plot.logp = function(minReads=4, minHetSamples=45) {
    
    opar = par(mfrow=c(10,1))
    par(mar=c(0.4,4,0.4,0.4))
    
    xmax = 305861025
    xmax = max(seg$pos[seg$contig==1])
    ymax = 8
    
    xlim = c(0,xmax)
    ylim = c(0,ymax)
    
    yline1 = rep(-log10(1e-2), 2)
    yline2 = rep(-log10(1e-6), 2)
    
    ## require a minimum number of Het samples
    filter = ((seg$Ahet+seg$Bhet)>=minHetSamples)

    ## require a minimum average number of reads per CALLED sample
    if (minReads>0) {
        filter = filter & ((seg$Aref+seg$Aalt+seg$Bref+seg$Balt)/(seg$Ahom+seg$Ahet+seg$Bhom+seg$Bhet)>=minReads)
    }

    ## total samples, hopefully
    nsamples = max(seg$Ahet+seg$Ahom+seg$Bhet+seg$Bhom)
    
    ## one plot per chromosome
    for (chr in 1:10) {
        
        plot(seg$pos[seg$contig==chr & filter], -log10(seg$pHetNC[seg$contig==chr & filter]), pch=1, col="black", ylab=paste("Chr",chr), xlim=xlim, ylim=ylim,  xaxt='n', xaxs='i')
        
        ## highlight significant p values
        points(seg$pos[seg$contig==chr & filter & seg$pHetNC<1e-2], -log10(seg$pHetNC[seg$contig==chr & filter & seg$pHetNC<1e-2]), pch=19, col="darkred")
        lines(xlim,  yline1, col="gray", lty=2)
        lines(xlim,  yline2, col="gray", lty=2)
        
        ## special cases
        if (chr==1) {
            ## title
            text(1, 0.85*ylim[2], cex=1.5, pos=4, "Fisher's Exact Test -log10(p) for Het/Ref vs Twin/Non-Twin")
        } else if (chr==2) {
            ## filter text
            subtitle = paste("Plotting filter: require at least",minHetSamples,"Het samples out of",nsamples)
            if (minReads>0) subtitle = paste(subtitle,"; loci must have at least",minReads,"average reads per called sample.")
            text(1, 0.85*ylim[2], cex=1.5, pos=4, subtitle)
        } else if (chr==3) {
            ## show ig1 location
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

}
