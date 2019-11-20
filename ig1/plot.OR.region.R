##
## plot odds ratio for each SNP on a given chromosomal region
##

plot.OR.region = function(chr, pos1, pos2, minReads=0, minHetSamples=45, log=FALSE, gaussFit=FALSE, fitRegion=c(0,0)) {

    ## red dots
    yhigh = 3
    
    ## some numbers
    nhet = seg$Ahet+seg$Bhet
    nreads = seg$Aref+seg$Aalt+seg$Bref+seg$Balt
    ncalled = seg$Ahom+seg$Ahet+seg$Bhom+seg$Bhet
    nsamples = max(ncalled) # hopefully the number of samples
    nAsamples = max(seg$Ahom+seg$Ahet) # hopefully the number of non-twin samples
    nBsamples = max(seg$Bhom+seg$Bhet) # hopefully the number of non-twin samples

    ## loci within given range
    filter = seg$contig==chr & seg$pos>=pos1 & seg$pos<=pos2

    ## require a minimum and maximum number of Het samples
    filter = filter & (nhet>=minHetSamples)

    ## require a minimum average number of reads per CALLED sample
    if (minReads>0) filter = filter & (nreads/ncalled>=minReads)

    ## Odds Ratio
    AhetOdds = seg$Ahet/(nAsamples-seg$Ahet)
    BhetOdds = seg$Bhet/(nBsamples-seg$Bhet)
    OR = BhetOdds[filter]/AhetOdds[filter]
    x = seg$pos[filter]

    ## THE PLOT
    ymax = max(OR[!is.nan(OR) & OR<1000])
    xlim = c(pos1,pos2)
    if (log) {
        plot(x, log10(OR), pch=1, col="black", xlab="Position", ylab="LOR", main=paste("Log Odds Ratio Het Twin / Het Non-Twin Chr",chr), xlim=xlim)
        lines(xlim, c(0,0), col="gray");
    } else {
        plot(x, OR, pch=1, col="black", xlab="Position", ylab="Odds Ratio", main=paste("Odds Ratio Het Twin / Het Non-Twin Chr",chr), xlim=xlim)
    }

    ## highlight high OR for twinning
    if (log) {
        points(x[OR>yhigh], log10(OR[OR>yhigh]), pch=19, col="darkred")
    } else {
        points(x[OR>yhigh], OR[OR>yhigh], pch=19, col="darkred")
    }

    ## title
    text(pos1, 0.95*ymax, cex=1.0, pos=4, "Odds Ratio for Het Twin / Het Non-Twin")
    ## subtitle
    subtitle = paste("Plotting filter: require at least",minHetSamples,"Het samples out of",nsamples)
    if (minReads>0) subtitle = paste(subtitle,"; loci must have at least",minReads,"average reads per called sample.")
    text(pos1, 0.90*ymax, cex=1.0, pos=4, subtitle)

    if (chr==3) {
        ## show ig1 location
        lines(c(171820386,171820386), c(0,ymax), col="black")
        lines(c(171823924,171823924), c(0,ymax), col="black")
        text(171822000, 0.75*ymax, "B73 LBD6 (ig1)", pos=4, col="black")
    } else if (chr==8) {
        ## marker umc1984
        points(c(81869662,81869383), c(0.75*ymax,0.75*ymax), col="darkorange", pch=2) 
        text(81869662, 0.75*ymax, "umc1984 (BE050190)", pos=4, col="darkorange")
    } else if (chr==9) {
        ## marker umc1586
        points(c(25250075,25250260), c(0.75*ymax,0.75*ymax), col="darkorange", pch=2) 
        text(25250075, 0.75*ymax, "umc1586 (AW053174)", pos=2, col="darkorange")
        ## marker umc2337
        points(c(26675959,26675253), c(0.75*ymax,0.75*ymax), col="darkorange", pch=2) 
        text(26675959, 0.75*ymax, "umc2337 (BM339026)", pos=4, col="darkorange")
    }

    ## fit a Gaussian to OR
    if (gaussFit) {
        if (max(fitRegion)>0) {
            y = OR[!is.nan(OR) & OR<1000 & OR>yhigh & x>=fitRegion[1] & x<=fitRegion[2]]
            x = x[!is.nan(OR) & OR<1000 & OR>yhigh & x>=fitRegion[1] & x<=fitRegion[2]]
        } else {
            y = OR[!is.nan(OR) & OR<1000 & OR>yhigh]
            x = x[!is.nan(OR) & OR<1000 & OR>yhigh]
        }
        mu = mean(x)
        sig = sd(x)
        scale = max(x)*max(y)
        fit = fitG(x, y, mu, sig, scale)
        yfit = fit$par[3]*dnorm(x,fit$par[1],fit$par[2])
        if (log) {
            lines(x, log10(yfit), col="blue", lwd=2)
        } else {
            lines(x, yfit, col="blue", lwd=2)
        }
        ks = ks.test(y,yfit)
        text(pos1, 0.85*ymax, cex=1.0, pos=4, paste("Gaussian fit peak (K-S p-value):",round(fit$par[1]),"(",ks["p.value"],")"), col="blue")
    }
    
}

fitG = function(x, y, mu, sig, scale) {
    f = function(p) {
        d = p[3]*dnorm(x,mean=p[1],sd=p[2])
        sum((d-y)^2)
    }
    optim(c(mu,sig,scale),f)
}

