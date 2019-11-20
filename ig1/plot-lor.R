##
## plot log odds ratio of each SNP on each chromosome, stacked plots
##

opar = par(mfrow=c(10,1))
par(mar=c(0.4,4,0.4,0.4))

maxLen = 305861025

xlim = c(0,maxLen)
ylim = c(-4,4)

## guide line
Aline = rep(-2,2)
Bline = rep( 2,2)


position = seg$position[seg$contig=="1"]
logOR = seg$logHetOR[seg$contig=="1"]
plot(position, logOR, ylab="Chr1", pch=1, xlim=xlim, ylim=ylim,  xaxt='n')
lines(xlim, c(0,0), col="black", lty=2)
lines(xlim, Aline, col="blue", lty=2)
lines(xlim, Bline, col="blue", lty=2)

position = seg$position[seg$contig=="2"]
logOR = seg$logHetOR[seg$contig=="2"]
plot(position, logOR, ylab="Chr2", pch=1, xlim=xlim, ylim=ylim,  xaxt='n')
lines(xlim, c(0,0), col="black", lty=2)
lines(xlim, Aline, col="blue", lty=2)
lines(xlim, Bline, col="blue", lty=2)

position = seg$position[seg$contig=="3"]
logOR = seg$logHetOR[seg$contig=="3"]
plot(position, logOR, ylab="Chr3", pch=1, xlim=xlim, ylim=ylim,  xaxt='n')
lines(xlim, c(0,0), col="black", lty=2)
lines(xlim, Aline, col="blue", lty=2)
lines(xlim, Bline, col="blue", lty=2)
lines(rep(180000000,2), ylim, col="red")
lines(rep(210000000,2), ylim, col="red")

position = seg$position[seg$contig=="4"]
logOR = seg$logHetOR[seg$contig=="4"]
plot(position, logOR, ylab="Chr4", pch=1, xlim=xlim, ylim=ylim,  xaxt='n')
lines(xlim, c(0,0), col="black", lty=2)
lines(xlim, Aline, col="blue", lty=2)
lines(xlim, Bline, col="blue", lty=2)

position = seg$position[seg$contig=="5"]
logOR = seg$logHetOR[seg$contig=="5"]
plot(position, logOR, ylab="Chr5", pch=1, xlim=xlim, ylim=ylim,  xaxt='n')
lines(xlim, c(0,0), col="black", lty=2)
lines(xlim, Aline, col="blue", lty=2)
lines(xlim, Bline, col="blue", lty=2)

position = seg$position[seg$contig=="6"]
logOR = seg$logHetOR[seg$contig=="6"]
plot(position, logOR, ylab="Chr6", pch=1, xlim=xlim, ylim=ylim,  xaxt='n')
lines(xlim, c(0,0), col="black", lty=2)
lines(xlim, Aline, col="blue", lty=2)
lines(xlim, Bline, col="blue", lty=2)

position = seg$position[seg$contig=="7"]
logOR = seg$logHetOR[seg$contig=="7"]
plot(position, logOR, ylab="Chr7", pch=1, xlim=xlim, ylim=ylim,  xaxt='n')
lines(xlim, c(0,0), col="black", lty=2)
lines(xlim, Aline, col="blue", lty=2)
lines(xlim, Bline, col="blue", lty=2)

position = seg$position[seg$contig=="8"]
logOR = seg$logHetOR[seg$contig=="8"]
plot(position, logOR, ylab="Chr8", pch=1, xlim=xlim, ylim=ylim,  xaxt='n')
lines(xlim, c(0,0), col="black", lty=2)
lines(xlim, Aline, col="blue", lty=2)
lines(xlim, Bline, col="blue", lty=2)

position = seg$position[seg$contig=="9"]
logOR = seg$logHetOR[seg$contig=="9"]
plot(position, logOR, ylab="Chr9", pch=1, xlim=xlim, ylim=ylim,  xaxt='n')
lines(xlim, c(0,0), col="black", lty=2)
lines(xlim, Aline, col="blue", lty=2)
lines(xlim, Bline, col="blue", lty=2)

position = seg$position[seg$contig=="10"]
logOR = seg$logHetOR[seg$contig=="10"]
plot(position, logOR, ylab="Chr10", pch=1, xlim=xlim, ylim=ylim,  xaxt='n')
lines(xlim, c(0,0), col="black", lty=2)
lines(xlim, Aline, col="blue", lty=2)
lines(xlim, Bline, col="blue", lty=2)

par(opar)
