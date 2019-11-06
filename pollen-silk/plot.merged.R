## dot plot the merged data frame.
##
## requires that load-merged.R be called in advance, which sets mixRatio for this set.

plot.merged = function(altAllele=FALSE) {

    if (altAllele) {
        allele = "ALT"
        mixRatio = altMixRatio
        deup = altdeup
        dedn = altdedn
        xData = merged$ALTtot.a
        yData = merged$ALTtot.m
        xSampleName = altSampleName
        ySampleName = mixSampleName
        mainTitle = paste(reference,"align.",mixSampleName,"ALT vs",altSampleName,"ALT @",refSampleName,"REF locs")
    } else {
        allele = "REF"
        mixRatio = refMixRatio
        deup = refdeup
        dedn = refdedn
        xData = merged$REFtot.r
        yData = merged$REFtot.m
        xSampleName = refSampleName
        ySampleName = mixSampleName
        mainTitle = paste(reference,"align.",mixSampleName,"REF vs",refSampleName,"REF @",altSampleName,"ALT locs")
    }

    ## DOT PLOT
    plot(xData[clean], yData[clean], 
         cex=0.5,
         log="xy",
         ## xlim=c(.5,100),
         ## ylim=c(100,20000),
         xlab = paste(xSampleName,allele,"count"),
         ylab = paste(ySampleName,allele,"count"),
         main=mainTitle
         )
    legend("topleft", bty="n",
           c(
               paste(singleTissues[xSampleName],"expression"),
               paste("Mix",allele,"ratio =",round(mixRatio,3))
           )
           )
           

    ## highlight up and dn
    points(xData[deup], yData[deup],
           cex=0.3, pch=19, col="darkgreen")
    points(xData[dedn], yData[dedn],
           cex=0.3, pch=19, col="darkred")

    ## print up and dn gene names
    text(xData[deup], yData[deup], merged$GENE[deup],
         pos=2, cex=0.5, offset=0.2, col="darkgreen")
    text(xData[dedn], yData[dedn], merged$GENE[dedn],
         pos=4, cex=0.5, offset=0.2, col="darkred")

    ## fixed ratio lines
    lines(c(1,1e5), c(1,1e5)*mixRatio, col="gray", lwd=3)
    lines(c(1,1e5), c(1,1e5)*mixRatio*2, col="darkgreen")
    lines(c(1,1e5), c(1,1e5)*mixRatio*4, col="darkgreen")
    lines(c(1,1e5), c(1,1e5)*mixRatio*8, col="darkgreen")
    lines(c(1,1e5), c(1,1e5)*mixRatio*16, col="darkgreen")
    lines(c(1,1e5), c(1,1e5)*mixRatio*32, col="darkgreen")
    lines(c(1,1e5), c(1,1e5)*mixRatio*64, col="darkgreen")
    lines(c(1,1e5), c(1,1e5)*mixRatio/2, col="darkred")
    lines(c(1,1e5), c(1,1e5)*mixRatio/4, col="darkred")
    lines(c(1,1e5), c(1,1e5)*mixRatio/8, col="darkred")
    lines(c(1,1e5), c(1,1e5)*mixRatio/16, col="darkred")
    lines(c(1,1e5), c(1,1e5)*mixRatio/32, col="darkred")
    lines(c(1,1e5), c(1,1e5)*mixRatio/64, col="darkred")

    ## ## print up and dn LOCI
    ## text(xData[deup], yData[deup],
    ##      paste(merged$CHR[deup],":",merged$POS[deup],sep=""),
    ##      pos=4, cex=0.7, offset=0.2, col="darkgreen")
    ## text(xData[dedn], yData[dedn],
    ##      paste(merged$CHR[dedn],":",merged$POS[dedn],sep=""),
    ##      pos=4, cex=0.7, offset=0.2, col="darkred")
}
