## dot plot the merged data frame.
## requires that mergVCF.R be called in advance.

## the scale factor to make MIX=REF for average loci.
# yscale = 3.94 
yscale = 0.549

quantiles = quantile(merged.clean$REFtot.m/merged.clean$REFtot.r, probs=c(0.01,0.1,0.5,0.9,0.99))
all = merged.clean$REFtot.m>0

highRatio = merged.clean$REFtot.m*yscale/merged.clean$REFtot.r > 16
lowRatio =  merged.clean$REFtot.m*yscale/merged.clean$REFtot.r < 1/16

selected = all

## DOT PLOT
plot(merged.clean$REFtot.r[selected], merged.clean$REFtot.m[selected], 
     cex=0.5,
     log="xy",
     xlim=c(1,1000),
     ylim=c(100,20000),
     xlab=paste(refSampleName,"(REF)"),
     ylab=paste(mixSampleName,"(MIX)"),
     main=paste("Alignment to",reference,mixSampleName,"vs",refSampleName,"@",altSampleName," ALT locs"))

## highlights
text(merged.clean$REFtot.r[highRatio], merged.clean$REFtot.m[highRatio],
     paste(merged.clean$CHR[highRatio],":",merged.clean$POS[highRatio],sep=""),
     pos=4, cex=0.5, offset=0.2, col="black")

text(merged.clean$REFtot.r[lowRatio], merged.clean$REFtot.m[lowRatio],
     paste(merged.clean$CHR[lowRatio],":",merged.clean$POS[lowRatio],sep=""),
     pos=4, cex=0.5, offset=0.2, col="black")

text(par()$xaxp[1], 10^(par()$usr[4]-0.2), paste("Lines at factors of two with mix scale factor =",yscale), pos=4, offset=0)

## points(merged.clean$REFtot.r[selected], merged.clean$REFtot.m[selected],
##        cex=0.4, pch=19, col="black")

## ratio lines
lines(c(1,1e5), c(1,1e5)/yscale, col="gray", lwd=3)
lines(c(1,1e5), c(1,1e5)/yscale*2, col="darkgreen")
lines(c(1,1e5), c(1,1e5)/yscale*4, col="darkgreen")
lines(c(1,1e5), c(1,1e5)/yscale*8, col="darkgreen")
lines(c(1,1e5), c(1,1e5)/yscale*16, col="darkgreen")
lines(c(1,1e5), c(1,1e5)/yscale*32, col="darkgreen")
lines(c(1,1e5), c(1,1e5)/yscale*64, col="darkgreen")
lines(c(1,1e5), c(1,1e5)/yscale/2, col="darkred")
lines(c(1,1e5), c(1,1e5)/yscale/4, col="darkred")
lines(c(1,1e5), c(1,1e5)/yscale/8, col="darkred")
lines(c(1,1e5), c(1,1e5)/yscale/16, col="darkred")
lines(c(1,1e5), c(1,1e5)/yscale/32, col="darkred")
lines(c(1,1e5), c(1,1e5)/yscale/64, col="darkred")

## ## RATIO PLOT
## plot(merged.clean$REFtot.r, merged.clean$REFtot.m/merged.clean$REFtot.r,
##      log="xy", cex=0.5,
##      xlim=c(1,3000), ylim=c(2,500),
##      xlab=paste(refSampleName,"(REF)"),
##      ylab=paste(mixSampleName,"(MIX)"),
##      main=paste("Alignment to",reference,":",mixSampleName, refSampleName," called @",altSampleName," ALT locs"))

## ## highlights
## text(merged.clean$REFtot.r[highRatio], merged.clean$REFtot.m[highRatio]/merged.clean$REFtot.r[highRatio],
##      paste(merged.clean$CHR[highRatio],":",merged.clean$POS[highRatio],sep=""), pos=4, cex=0.5, offset=0.2, col="black")
## points(merged.clean$REFtot.r[highRatio], merged.clean$REFtot.m[highRatio]/merged.clean$REFtot.r[highRatio],
##        col="darkgreen", cex=0.4, pch=19)

## text(merged.clean$REFtot.r[lowRatio], merged.clean$REFtot.m[lowRatio]/merged.clean$REFtot.r[lowRatio],
##      paste(merged.clean$CHR[lowRatio],":",merged.clean$POS[lowRatio],sep=""), pos=4, cex=0.5, offset=0.2, col="black")
## points(merged.clean$REFtot.r[lowRatio], merged.clean$REFtot.m[lowRatio]/merged.clean$REFtot.r[lowRatio],
##        col="darkred", cex=0.4, pch=19)

## lines(c(1,1e5), c(1,1)*quantiles[[1]], col="gray")
## lines(c(1,1e5), c(1,1)*quantiles[[2]], col="gray", lwd=2)


