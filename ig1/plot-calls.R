##
## plot the number of genotype calls for each group per chromosome
##

opar = par(mfrow=c(10,1))
par(mar=c(0.4,4,0.4,0.4))

maxLen = 305861025

xlim = c(0,maxLen)
ylim = c(-50,50)

yline = c(24,24)

plot(seg$position[seg$contig=="1"], seg$Bhet[seg$contig=="1"],  ylab="Chr1", pch=1, cex=0.5, col="darkgreen", xlim=xlim, ylim=ylim,  xaxt='n', xaxs="i",yaxs="i")
points(seg$position[seg$contig=="1"], -(seg$Bhom[seg$contig=="1"]), pch=1, cex=0.5, col="darkblue")
lines(xlim, yline, col="red")
lines(xlim,-yline, col="red")

text(0, -40, "Twinning Group", pos=4, cex=2)

plot(seg$position[seg$contig=="2"], seg$Bhet[seg$contig=="2"], ylab="Chr2", pch=1, cex=0.5, col="darkgreen", xlim=xlim, ylim=ylim,  xaxt='n', xaxs="i",yaxs="i")
points(seg$position[seg$contig=="2"], -(seg$Bhom[seg$contig=="2"]), pch=1, cex=0.5, col="darkblue")
lines(xlim, yline, col="red")
lines(xlim,-yline, col="red")

plot(seg$position[seg$contig=="3"], seg$Bhet[seg$contig=="3"], ylab="Chr3", pch=1, cex=0.5, col="darkgreen", xlim=xlim, ylim=ylim,  xaxt='n', xaxs="i",yaxs="i")
points(seg$position[seg$contig=="3"], -(seg$Bhom[seg$contig=="3"]), pch=1, cex=0.5, col="darkblue")
lines(xlim, yline, col="red")
lines(xlim,-yline, col="red")

plot(seg$position[seg$contig=="4"], seg$Bhet[seg$contig=="4"], ylab="Chr4", pch=1, cex=0.5, col="darkgreen", xlim=xlim, ylim=ylim,  xaxt='n', xaxs="i",yaxs="i")
points(seg$position[seg$contig=="4"], -(seg$Bhom[seg$contig=="4"]), pch=1, cex=0.5, col="darkblue")
lines(xlim, yline, col="red")
lines(xlim,-yline, col="red")

plot(seg$position[seg$contig=="5"], seg$Bhet[seg$contig=="5"], ylab="Chr5", pch=1, cex=0.5, col="darkgreen", xlim=xlim, ylim=ylim,  xaxt='n', xaxs="i",yaxs="i")
points(seg$position[seg$contig=="5"], -(seg$Bhom[seg$contig=="5"]), pch=1, cex=0.5, col="darkblue")
lines(xlim, yline, col="red")
lines(xlim,-yline, col="red")

plot(seg$position[seg$contig=="6"], seg$Bhet[seg$contig=="6"], ylab="Chr6", pch=1, cex=0.5, col="darkgreen", xlim=xlim, ylim=ylim,  xaxt='n', xaxs="i",yaxs="i")
points(seg$position[seg$contig=="6"], -(seg$Bhom[seg$contig=="6"]), pch=1, cex=0.5, col="darkblue")
lines(xlim, yline, col="red")
lines(xlim,-yline, col="red")

plot(seg$position[seg$contig=="7"], seg$Bhet[seg$contig=="7"], ylab="Chr7", pch=1, cex=0.5, col="darkgreen", xlim=xlim, ylim=ylim,  xaxt='n', xaxs="i",yaxs="i")
points(seg$position[seg$contig=="7"], -(seg$Bhom[seg$contig=="7"]), pch=1, cex=0.5, col="darkblue")
lines(xlim, yline, col="red")
lines(xlim,-yline, col="red")

plot(seg$position[seg$contig=="8"], seg$Bhet[seg$contig=="8"], ylab="Chr8", pch=1, cex=0.5, col="darkgreen", xlim=xlim, ylim=ylim,  xaxt='n', xaxs="i",yaxs="i")
points(seg$position[seg$contig=="8"], -(seg$Bhom[seg$contig=="8"]), pch=1, cex=0.5, col="darkblue")
lines(xlim, yline, col="red")
lines(xlim,-yline, col="red")

plot(seg$position[seg$contig=="9"], seg$Bhet[seg$contig=="9"], ylab="Chr9", pch=1, cex=0.5, col="darkgreen", xlim=xlim, ylim=ylim,  xaxt='n', xaxs="i",yaxs="i")
points(seg$position[seg$contig=="9"], -(seg$Bhom[seg$contig=="9"]), pch=1, cex=0.5, col="darkblue")
lines(xlim, yline, col="red")
lines(xlim,-yline, col="red")

plot(seg$position[seg$contig=="10"], seg$Bhet[seg$contig=="10"], ylab="Chr10", pch=1, cex=0.5, col="darkgreen", xlim=xlim, ylim=ylim,  xaxt='n', xaxs="i",yaxs="i")
points(seg$position[seg$contig=="10"],-(+seg$Bhom[seg$contig=="10"]), pch=1, cex=0.5, col="darkblue")
lines(xlim, yline, col="red")
lines(xlim,-yline, col="red")

par(opar)
