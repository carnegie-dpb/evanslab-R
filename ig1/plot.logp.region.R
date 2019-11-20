##
## plot log Fisher p for each SNP on a given chromosomal region
##

plot.logp.region = function(chr, pos1, pos2, minReads=0, minHetSamples=45, gaussFit=FALSE, fitRegion=c(0,0)) {

    ymax = 8
    
    xlim = c(pos1,pos2)
    ylim = c(0,ymax)
    
    yline1 = rep(-log10(1e-2), 2)
    yline2 = rep(-log10(1e-6), 2)

    ## require a minimum and maximum number of Het samples
    filter = ((seg$Ahet+seg$Bhet)>=minHetSamples)

    ## require a minimum average number of reads per CALLED sample
    if (minReads>0) {
        filter = filter & ((seg$Aref+seg$Aalt+seg$Bref+seg$Balt)/(seg$Ahom+seg$Ahet+seg$Bhom+seg$Bhet)>=minReads)
    }

    ## within given range
    filter = filter & seg$contig==chr & seg$pos>=pos1 & seg$pos<=pos2

    ## total samples, hopefully
    nsamples = max(seg$Ahet+seg$Ahom+seg$Bhet+seg$Bhom)

    ## PLOT
    plot(seg$pos[filter], -log10(seg$pHetNC[filter]), pch=1, col="black", xlab="Position", ylab="-log10(p)", main=paste("Fisher's Exact Test Chr",chr), xlim=xlim, ylim=ylim)

    ## highlight significant p values
    points(seg$pos[filter & seg$pHetNC<1e-2], -log10(seg$pHetNC[filter & seg$pHetNC<1e-2]), pch=19, col="darkred")

    ## significance lines
    lines(xlim,  yline1, col="gray", lty=2)
    lines(xlim,  yline2, col="gray", lty=2)

    ## title
    text(pos1, 0.95*ylim[2], cex=1.0, pos=4, "Fisher's Exact Test -log10(p) for Het/Ref vs Twin/Non-Twin")
    ## subtitle
    subtitle = paste("Plotting filter: require at least",minHetSamples,"Het samples out of",nsamples)
    if (minReads>0) subtitle = paste(subtitle,"; loci must have at least",minReads,"average reads per called sample.")
    text(pos1, 0.90*ylim[2], cex=1.0, pos=4, subtitle)

    if (chr==1) {
        ## OsO3g0241300 BLAST locations
        lines(c(32818769,32819151), 0.8*c(0,ymax), col="black")
        lines(c(32816916,32817121), 0.8*c(0,ymax), col="black")
        text(32818769, 0.8*ymax, "OsPE BLAST match", col="black")
    } else if (chr==3) {
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
        ## OsO3g0241300 BLAST locations
        lines(c(146099175,146099584), 0.8*c(0,ymax), col="black")
        lines(c(146100542,146100746), 0.8*c(0,ymax), col="black")
        text(146099175, 0.8*ymax, "OsPE BLAST match", col="black")
    }

    ## fit a Gaussian
    if (gaussFit) {
        if (max(fitRegion)>0) {
            x = seg$pos[filter & seg$pHetNC<1e-2 & seg$pos>=fitRegion[1] & seg$pos<=fitRegion[2]]
            y = -log10(seg$pHetNC[filter & seg$pHetNC<1e-2 & seg$pos>=fitRegion[1] & seg$pos<=fitRegion[2]])
        } else {
            x = seg$pos[filter & seg$pHetNC<1e-2]
            y = -log10(seg$pHetNC[filter & seg$pHetNC<1e-2])
        }
        mu = mean(x)
        sig = sd(x)
        scale = max(x)*max(y)
        fit = fitG(x, y, mu, sig, scale)
        yfit = fit$par[3]*dnorm(x,fit$par[1],fit$par[2])
        lines(x, yfit, col="blue", lwd=2)
        ks = ks.test(y,yfit)
        text(pos1, 0.85*ylim[2], cex=1.0, pos=4, paste("Gaussian fit peak (K-S p-value):",round(fit$par[1]),"(",ks["p.value"],")"), col="blue")
    }
    
}

fitG = function(x, y, mu, sig, scale) {
    f = function(p) {
        d = p[3]*dnorm(x,mean=p[1],sd=p[2])
        sum((d-y)^2)
    }
    optim(c(mu,sig,scale),f)
}
