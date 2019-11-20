##
## plot log Fisher p for each SNP on a given pair of chromosomal regions, and mark syntenic regions from the dataframe
##
## df: length idents gaps bits eval chrAstart chrAend chrBstart chrBend chrA chrB

plot.logp.tworegions = function(df, chrA, chrAstart, chrAend, chrB, chrBstart, chrBend, minReads=0, minHetSamples=45) {

    ymax = 8

    xlimA = c(chrAstart, chrAend)
    xlimB = c(chrBstart, chrBend)

    ylim = c(0,ymax)
    
    yline1 = rep(-log10(1e-2), 2)
    yline2 = rep(-log10(1e-6), 2)

    ## require a minimum and maximum number of Het samples
    filter = ((seg$Ahet+seg$Bhet)>=minHetSamples)

    ## require a minimum average number of reads per CALLED sample
    if (minReads>0) {
        filter = filter & ((seg$Aref+seg$Aalt+seg$Bref+seg$Balt)/(seg$Ahom+seg$Ahet+seg$Bhom+seg$Bhet)>=minReads)
    }

    ## within given range on chrA, chrB
    filterA = filter & seg$contig==chrA & seg$pos>=chrAstart & seg$pos<=chrAend
    filterB = filter & seg$contig==chrB & seg$pos>=chrBstart & seg$pos<=chrBend

    ## total samples, hopefully
    nsamples = max(seg$Ahet+seg$Ahom+seg$Bhet+seg$Bhom)

    options(scipen=5)
    oldpar = par()
    par(
        mfrow=(c(2,1)),
        mar=c(2,4,1,1)+0.1 # bottom,left,top,right
    )

    ## PLOT A
    plot(seg$pos[filterA], -log10(seg$pHetNC[filterA]), pch=1, col="darkgreen", xlab="", ylab=paste("Chr",chrA,"-log10(p)"), xlim=xlimA, ylim=ylim, cex=0.5)
    points(seg$pos[filterA & seg$pHetNC<1e-2], -log10(seg$pHetNC[filterA & seg$pHetNC<1e-2]), pch=19, col="darkgreen", cex=0.5)
    lines(xlimA,  yline1, col="gray", lty=2)
    lines(xlimA,  yline2, col="gray", lty=2)
    for (i in 1:50) {
        lines(c(df$chrAstart[i],df$chrAend[i]), c(ymax,ymax)*0.90, col="black")
        text((df$chrAstart[i]+df$chrAend[i])/2, ymax*0.95, i, col="black", cex=0.75)
    }

    plot(seg$pos[filterB], -log10(seg$pHetNC[filterB]), pch=1, col="darkblue",  xlab="", ylab=paste("Chr",chrB,"-log10(p)"), xlim=xlimB, ylim=ylim, cex=0.5)
    points(seg$pos[filterB & seg$pHetNC<1e-2], -log10(seg$pHetNC[filterB & seg$pHetNC<1e-2]), pch=19, col="darkblue", cex=0.5)
    lines(xlimB,  yline1, col="gray", lty=2)
    lines(xlimB,  yline2, col="gray", lty=2)
    for (i in 1:50) {
        lines(c(df$chrBstart[i],df$chrBend[i]), c(ymax,ymax)*0.90, col="black")
        text((df$chrBstart[i]+df$chrBend[i])/2, ymax*0.95, i, col="black", cex=0.75)
    }


    ##     ## title
    ##     text(pos1, 0.95*ylim[2], cex=1.0, pos=4, "Fisher's Exact Test -log10(p) for Het/Ref vs Twin/Non-Twin")
    ##     ## subtitle
    ##     subtitle = paste("Plotting filter: require at least",minHetSamples,"Het samples out of",nsamples)
    ##     if (minReads>0) subtitle = paste(subtitle,"; loci must have at least",minReads,"average reads per called sample.")
    ##     text(pos1, 0.90*ylim[2], cex=1.0, pos=4, subtitle)

    par = oldpar
    
}

