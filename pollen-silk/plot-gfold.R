## plot the GFOLD diff results
## GeneSymbol	GeneName	GFOLD(0.01)	E-FDR	log2fdc	1stRPKM	2ndRPKM

plot(gfold$RPKM1/norm1, gfold$RPKM2/norm2,
     xlab=paste(sample1,"(RPKM)"), ylab=paste(sample2,"(norm. RPKM)"),
     main=paste(reference,"perfect matches:",sample2,"vs",sample1),
     ## xlim=c(0,1),
     ## ylim=c(0,1),
     log="xy",
     cex=0.5
     )

gfoldUp = gfold$GFOLD >  4.0
gfoldDn = gfold$GFOLD < -4.0

points(gfold$RPKM1[gfoldUp]/norm1, gfold$RPKM2[gfoldUp]/norm2,
       pch=19, cex=0.3, col="darkgreen"
       )

points(gfold$RPKM1[gfoldDn]/norm1, gfold$RPKM2[gfoldDn]/norm2,
       pch=19, cex=0.3, col="darkred"
       )

lines(c(1e-3,1e5), c(1e-3,1e5), col="gray", lwd=3)
lines(c(1e-3,1e5), c(1e-3,1e5)*2, col="darkgreen", lwd=1)
lines(c(1e-3,1e5), c(1e-3,1e5)*4, col="darkgreen", lwd=1)
lines(c(1e-3,1e5), c(1e-3,1e5)*8, col="darkgreen", lwd=1)
lines(c(1e-3,1e5), c(1e-3,1e5)*16, col="darkgreen", lwd=1)
lines(c(1e-3,1e5), c(1e-3,1e5)/2, col="darkred", lwd=1)
lines(c(1e-3,1e5), c(1e-3,1e5)/4, col="darkred", lwd=1)
lines(c(1e-3,1e5), c(1e-3,1e5)/8, col="darkred", lwd=1)
lines(c(1e-3,1e5), c(1e-3,1e5)/16, col="darkred", lwd=1)

text(par()$xaxp[1], par()$yaxp[2], paste("Lines at factors of two.","\n",sample2,"/",sample1,"=",round(norm2,3), sep=""), adj=c(0,1))

text(gfold$RPKM1[gfoldUp], gfold$RPKM2[gfoldUp], rownames(gfold)[gfoldUp], 
     cex=0.6, pos=4, offset=0.2, col="darkgreen")
text(gfold$RPKM1[gfoldDn], gfold$RPKM2[gfoldDn], rownames(gfold)[gfoldDn], 
     cex=0.6, pos=4, offset=0.2, col="darkred")
