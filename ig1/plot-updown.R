##
## plot sample counts for each SNP, up for B group, down for A group
##

## opar = par(mfrow=c(10,1))
## par(mar=c(0.4,4,0.4,0.4))

Atot = 45
Btot = 48

maxLen = 305861025

xlim = c(0,maxLen)

plot(seg.both$position[seg.both$contig=="1"],    seg.both$Bcount[seg.both$contig=="1"], ylab="Chr1", pch=19, xlim=xlim, xaxt='n', ylim=c(-45,48), col="red")
points(seg.both$position[seg.both$contig=="1"], -seg.both$Acount[seg.both$contig=="1"], pch=19, col="blue")
lines(xlim, rep(  3,2), col="black")

## plot(seg.both$position[seg.both$contig=="2"], seg.both$Bcount[seg.both$contig=="2"]-seg.both$Acount[seg.both$contig=="2"], ylab="Chr2", pch=19, xlim=xlim, xaxt='n', ylim=c(-30,30))
## lines(xlim, rep(  3,2), col="black")
## lines(xlim, rep( 23,2), col="red")
## lines(xlim, rep(-17,2), col="red")

## plot(seg.both$position[seg.both$contig=="3"], seg.both$Bcount[seg.both$contig=="3"]-seg.both$Acount[seg.both$contig=="3"], ylab="Chr3", pch=19, xlim=xlim, xaxt='n', ylim=c(-30,30))
## lines(xlim, rep(  3,2), col="black")
## lines(xlim, rep( 23,2), col="red")
## lines(xlim, rep(-17,2), col="red")

## plot(seg.both$position[seg.both$contig=="4"], seg.both$Bcount[seg.both$contig=="4"]-seg.both$Acount[seg.both$contig=="4"], ylab="Chr4", pch=19, xlim=xlim, xaxt='n', ylim=c(-30,30))
## lines(xlim, rep(  3,2), col="black")
## lines(xlim, rep( 23,2), col="red")
## lines(xlim, rep(-17,2), col="red")

## plot(seg.both$position[seg.both$contig=="5"], seg.both$Bcount[seg.both$contig=="5"]-seg.both$Acount[seg.both$contig=="5"], ylab="Chr5", pch=19, xlim=xlim, xaxt='n', ylim=c(-30,30))
## lines(xlim, rep(  3,2), col="black")
## lines(xlim, rep( 23,2), col="red")
## lines(xlim, rep(-17,2), col="red")

## plot(seg.both$position[seg.both$contig=="6"], seg.both$Bcount[seg.both$contig=="6"]-seg.both$Acount[seg.both$contig=="6"], ylab="Chr6", pch=19, xlim=xlim, xaxt='n', ylim=c(-30,30))
## lines(xlim, rep(  3,2), col="black")
## lines(xlim, rep( 23,2), col="red")
## lines(xlim, rep(-17,2), col="red")

## plot(seg.both$position[seg.both$contig=="7"], seg.both$Bcount[seg.both$contig=="7"]-seg.both$Acount[seg.both$contig=="7"], ylab="Chr7", pch=19, xlim=xlim, xaxt='n', ylim=c(-30,30))
## lines(xlim, rep(  3,2), col="black")
## lines(xlim, rep( 23,2), col="red")
## lines(xlim, rep(-17,2), col="red")

## plot(seg.both$position[seg.both$contig=="8"], seg.both$Bcount[seg.both$contig=="8"]-seg.both$Acount[seg.both$contig=="8"], ylab="Chr8", pch=19, xlim=xlim, xaxt='n', ylim=c(-30,30))
## lines(xlim, rep(  3,2), col="black")
## lines(xlim, rep( 23,2), col="red")
## lines(xlim, rep(-17,2), col="red")

## plot(seg.both$position[seg.both$contig=="9"], seg.both$Bcount[seg.both$contig=="9"]-seg.both$Acount[seg.both$contig=="9"], ylab="Chr9", pch=19, xlim=xlim, xaxt='n', ylim=c(-30,30))
## lines(xlim, rep(  3,2), col="black")
## lines(xlim, rep( 23,2), col="red")
## lines(xlim, rep(-17,2), col="red")

## plot(seg.both$position[seg.both$contig=="10"], seg.both$Bcount[seg.both$contig=="10"]-seg.both$Acount[seg.both$contig=="10"], ylab="Chr10", pch=19, xlim=xlim, xaxt='n', ylim=c(-30,30))
## lines(xlim, rep(  3,2), col="black")
## lines(xlim, rep( 23,2), col="red")
## lines(xlim, rep(-17,2), col="red")

## par(opar)
