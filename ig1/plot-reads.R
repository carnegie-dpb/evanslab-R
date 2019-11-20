##
## plot the ref and alt reads for a given group
##

opar = par(mfrow=c(10,1))
par(mar=c(0.4,4,0.4,0.4))

maxLen = 305861025

xlim = c(0,maxLen)
ylim = c(-1,1)

yline = c(24,24)

plot(seg$position[seg$contig=="1"], (seg$Balt[seg$contig=="1"]-seg$Bref[seg$contig=="1"])/(seg$Balt[seg$contig=="1"]+seg$Bref[seg$contig=="1"]),
     ylab="Chr1", pch=1, cex=0.5, col="black", xlim=xlim, ylim=ylim,  xaxt='n', xaxs="i",yaxs="i")
lines(xlim, rep(1,2), col="red", lty=2)
lines(xlim, rep(0,2), col="red", lty=2)

text(0, 0, "Twinning Group", pos=4, cex=2, col="white")

plot(seg$position[seg$contig=="2"], (seg$Balt[seg$contig=="2"]-seg$Bref[seg$contig=="2"])/(seg$Balt[seg$contig=="2"]+seg$Bref[seg$contig=="2"]),
     ylab="Chr2", pch=1, cex=0.5, col="black", xlim=xlim, ylim=ylim,  xaxt='n', xaxs="i",yaxs="i")
lines(xlim, rep(1,2), col="red", lty=2)
lines(xlim, rep(0,2), col="red", lty=2)

plot(seg$position[seg$contig=="3"], (seg$Balt[seg$contig=="3"]-seg$Bref[seg$contig=="3"])/(seg$Balt[seg$contig=="3"]+seg$Bref[seg$contig=="3"]),
     ylab="Chr3", pch=1, cex=0.5, col="black", xlim=xlim, ylim=ylim,  xaxt='n', xaxs="i",yaxs="i")
lines(xlim, rep(1,2), col="red", lty=2)
lines(xlim, rep(0,2), col="red", lty=2)

plot(seg$position[seg$contig=="4"], (seg$Balt[seg$contig=="4"]-seg$Bref[seg$contig=="4"])/(seg$Balt[seg$contig=="4"]+seg$Bref[seg$contig=="4"]),
     ylab="Chr4", pch=1, cex=0.5, col="black", xlim=xlim, ylim=ylim,  xaxt='n', xaxs="i",yaxs="i")
lines(xlim, rep(1,2), col="red", lty=2)
lines(xlim, rep(0,2), col="red", lty=2)

plot(seg$position[seg$contig=="5"], (seg$Balt[seg$contig=="5"]-seg$Bref[seg$contig=="5"])/(seg$Balt[seg$contig=="5"]+seg$Bref[seg$contig=="5"]),
     ylab="Chr5", pch=1, cex=0.5, col="black", xlim=xlim, ylim=ylim,  xaxt='n', xaxs="i",yaxs="i")
lines(xlim, rep(1,2), col="red", lty=2)
lines(xlim, rep(0,2), col="red", lty=2)

plot(seg$position[seg$contig=="6"], (seg$Balt[seg$contig=="6"]-seg$Bref[seg$contig=="6"])/(seg$Balt[seg$contig=="6"]+seg$Bref[seg$contig=="6"]),
     ylab="Chr6", pch=1, cex=0.5, col="black", xlim=xlim, ylim=ylim,  xaxt='n', xaxs="i",yaxs="i")
lines(xlim, rep(1,2), col="red", lty=2)
lines(xlim, rep(0,2), col="red", lty=2)

plot(seg$position[seg$contig=="7"], (seg$Balt[seg$contig=="7"]-seg$Bref[seg$contig=="7"])/(seg$Balt[seg$contig=="7"]+seg$Bref[seg$contig=="7"]),
     ylab="Chr7", pch=1, cex=0.5, col="black", xlim=xlim, ylim=ylim,  xaxt='n', xaxs="i",yaxs="i")
lines(xlim, rep(1,2), col="red", lty=2)
lines(xlim, rep(0,2), col="red", lty=2)

plot(seg$position[seg$contig=="8"], (seg$Balt[seg$contig=="8"]-seg$Bref[seg$contig=="8"])/(seg$Balt[seg$contig=="8"]+seg$Bref[seg$contig=="8"]),
     ylab="Chr8", pch=1, cex=0.5, col="black", xlim=xlim, ylim=ylim,  xaxt='n', xaxs="i",yaxs="i")
lines(xlim, rep(1,2), col="red", lty=2)
lines(xlim, rep(0,2), col="red", lty=2)

plot(seg$position[seg$contig=="9"], (seg$Balt[seg$contig=="9"]-seg$Bref[seg$contig=="9"])/(seg$Balt[seg$contig=="9"]+seg$Bref[seg$contig=="9"]),
     ylab="Chr9", pch=1, cex=0.5, col="black", xlim=xlim, ylim=ylim,  xaxt='n', xaxs="i",yaxs="i")
lines(xlim, rep(1,2), col="red", lty=2)
lines(xlim, rep(0,2), col="red", lty=2)

plot(seg$position[seg$contig=="10"], (seg$Balt[seg$contig=="10"]-seg$Bref[seg$contig=="10"])/(seg$Balt[seg$contig=="10"]+seg$Bref[seg$contig=="10"]),
     ylab="Chr10", pch=1, cex=0.5, col="black", xlim=xlim, ylim=ylim,  xaxt='n', xaxs="i",yaxs="i")
lines(xlim, rep(1,2), col="red", lty=2)
lines(xlim, rep(0,2), col="red", lty=2)


par(opar)
