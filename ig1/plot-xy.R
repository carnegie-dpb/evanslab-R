##
## x-y plot the count or reads of A vs B with guiding lines
##

plot(seg.both$Acount[seg.both$contig=="1"], seg.both$Bcount[seg.both$contig=="1"], xlab="A count", ylab="B count", pch=19, xlim=c(0,48), ylim=c(0,48))
lines(c(0,48),c(0,48), col="red")
lines(c(0,48),c(0,48)*2, col="red")
lines(c(0,48),c(0,48)*4, col="red")
lines(c(0,48),c(0,48)*8, col="red")
lines(c(0,48),c(0,48)*16, col="red")
