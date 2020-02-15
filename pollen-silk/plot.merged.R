## dot plot the merged data frame.
##
## requires that read-merged.R or load-merged.R be called in advance, which sets mixRatio for this set.

plot.merged = function(altAllele=FALSE, upGene=NULL, dnGene=NULL, highlight=FALSE, xmin=3, xmax=1e4, ymin=3, ymax=1e4) {

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
        mainTitle = paste(reference, "alignment:",
                          sampleTissues[mixSampleName],"REF vs",
                          sampleTissues[refSampleName],"REF called @",
                          sampleTissues[altSampleName],"ALT locs")
    }

    ## DOT PLOT
    plot(xData[clean], yData[clean], 
         cex=0.5,
         log="xy",
         xlim=c(xmin, xmax),
         ylim=c(ymin, ymax),
         xlab = paste(sampleTissues[xSampleName],allele,"count"),
         ylab = paste(sampleTissues[ySampleName],allele,"count")
         ## main=mainTitle, cex.main=0.75
         )

    ## fixed ratio lines
    lines(c(1,1e5), c(1,1e5)*mixRatio, col="gray", lwd=3)
    lines(c(1,1e5), c(1,1e5)*mixRatio*2, col="darkred")
    lines(c(1,1e5), c(1,1e5)*mixRatio*4, col="darkred")
    lines(c(1,1e5), c(1,1e5)*mixRatio*8, col="darkred")
    lines(c(1,1e5), c(1,1e5)*mixRatio*16, col="darkred")
    lines(c(1,1e5), c(1,1e5)*mixRatio*32, col="darkred")
    lines(c(1,1e5), c(1,1e5)*mixRatio*64, col="darkred")
    lines(c(1,1e5), c(1,1e5)*mixRatio/2, col="darkblue")
    lines(c(1,1e5), c(1,1e5)*mixRatio/4, col="darkblue")
    lines(c(1,1e5), c(1,1e5)*mixRatio/8, col="darkblue")
    lines(c(1,1e5), c(1,1e5)*mixRatio/16, col="darkblue")
    lines(c(1,1e5), c(1,1e5)*mixRatio/32, col="darkblue")
    lines(c(1,1e5), c(1,1e5)*mixRatio/64, col="darkblue")

    ## show points for a given UP gene
    if (!is.null(upGene)) {
        genepts = merged$GENE==upGene
        points(xData[genepts], yData[genepts], cex=1.2, pch=19, col="darkred")
        if (reference=="W22") {
            synpts = B73.W22.syntelogs$W22_gene_model_id==upGene
            syntelogs = B73.W22.syntelogs$B73v4_gene_model_id[synpts]
        } else {
            synpts = B73.W22.syntelogs$B73v4_gene_model_id==upGene
            syntelogs = B73.W22.syntelogs$W22_gene_model_id[synpts]
        }
    }

    ## show points for a given DN gene
    if (!is.null(dnGene)) {
        genepts = merged$GENE==dnGene
        points(xData[genepts], yData[genepts], cex=1.2, pch=19, col="darkblue")
        if (reference=="W22") {
            synpts = B73.W22.syntelogs$W22_gene_model_id==dnGene
            syntelogs = B73.W22.syntelogs$B73v4_gene_model_id[synpts]
        } else {
            synpts = B73.W22.syntelogs$B73v4_gene_model_id==dnGene
            syntelogs = B73.W22.syntelogs$W22_gene_model_id[synpts]
        }
    }

    text(xmax, ymin, paste("mix/single",allele,"ratio =",round(mixRatio,3)), pos=2)
    
    ## legend in top left
    legend = c("no change",
               "×2, ×4, ×8, ...",
               "×1/2, ×1/4, ×1/8, ..."
               )
    if (!is.null(upGene)) {
        legend = c(legend, upGene)
    }
    if (!is.null(dnGene)) {
        legend = c(legend, dnGene)
    }
    
    legend("topleft", bty="n",
           legend,
           pch=c(-1,-1,-1,19,19),
           lty=c(1,1,1,0,0),
           lwd=c(3,1,1,0,0),
           pt.cex=c(0,0,0,1.2,1.2),
           col=c("darkgray","darkred","darkblue","darkred","darkblue"),
           text.col=c("darkgray","darkred","darkblue","darkred","darkblue")
           )

    ## highlight up and dn
    if (highlight) {
        points(xData[deup], yData[deup],
               cex=0.3, pch=19, col="darkred")
        points(xData[dedn], yData[dedn],
               cex=0.3, pch=19, col="darkblue")

        ## print up and dn gene names
        text(xData[deup], yData[deup], merged$GENE[deup],
             pos=4, cex=0.6, offset=0.2, col="darkred")
        text(xData[dedn], yData[dedn], merged$GENE[dedn],
             pos=4, cex=0.6, offset=0.2, col="darkblue")
    }

    ## ## print up and dn LOCI
    ## text(xData[deup], yData[deup],
    ##      paste(merged$CHR[deup],":",merged$POS[deup],sep=""),
    ##      pos=4, cex=0.7, offset=0.2, col="darkred")
    ## text(xData[dedn], yData[dedn],
    ##      paste(merged$CHR[dedn],":",merged$POS[dedn],sep=""),
    ##      pos=4, cex=0.7, offset=0.2, col="darkblue")
}
