##
## do a dotplot to divine loci that are significantly differently expressed between two samples
##
yscale = 1.0
dotplot = function(xcounts, ycounts) {
    min = 1e1
    max = 1e6
    xmed = median(xcounts)
    ymed = median(ycounts)
    xnorm = (xmed+ymed)/xmed/2
    ynorm = (xmed+ymed)/ymed/2
    plot(xcounts, ycounts*yscale,
         xlab=deparse(substitute(xcounts)), ylab=deparse(substitute(ycounts)),
         xlim=c(min,max), ylim=c(min,max),
         cex=0.5, log="xy", col="black")
    text(xcounts, ycounts, rownames(htseq.B73), cex=0.5, pos=4, col="black")
    lines(c(min,max), c(min,max), col="red", lwd=1)
}
