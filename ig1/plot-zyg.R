##
## plot log odds ratio of each group's hom/het SNP on each chromosome, stacked plots
##

opar = par(mfrow=c(10,1))
par(mar=c(0.4,4,0.4,0.4))

Atot = 45
Btot = 48

maxLen = 305861025

xlim = c(0,maxLen)
ylim = c(-10,10)

## lines at 20 Hom:0 Het and vice versa
Nline = 20
Aline = log( (Nline+0.5)/(A-Nline+0.5) / ( 0.5/(A+0.5) ) )
Bline = log( (Nline+0.5)/(B-Nline+0.5) / ( 0.5/(B+0.5) ) )

## contig position ref alt Areads Breads Ahom Ahet   lnAOR Bhom Bhet   lnBOR

plot(zyg$position[zyg$contig=="1"], zyg$lnBOR[zyg$contig=="1"], ylab="Chr1", pch=1, cex=0.5, col="red", xlim=xlim, ylim=ylim,  xaxt='n', xaxs="i",yaxs="i")
points(zyg$position[zyg$contig=="1"], zyg$lnAOR[zyg$contig=="1"], pch=1, cex=0.5, col="green")
lines(xlim, rep( Aline,2), lty=3, col="red")
lines(xlim, rep(-Aline,2), lty=3, col="red")
    
plot(zyg$position[zyg$contig=="2"], zyg$lnBOR[zyg$contig=="2"], ylab="Chr2", pch=1, cex=0.5, col="red", xlim=xlim, ylim=ylim,  xaxt='n', xaxs="i",yaxs="i")
points(zyg$position[zyg$contig=="2"], zyg$lnAOR[zyg$contig=="2"], pch=1, cex=0.5, col="green")
lines(xlim, rep( Aline,2), lty=3, col="red")
lines(xlim, rep(-Aline,2), lty=3, col="red")

plot(zyg$position[zyg$contig=="3"], zyg$lnBOR[zyg$contig=="3"], ylab="Chr3", pch=1, cex=0.5, col="red", xlim=xlim, ylim=ylim,  xaxt='n', xaxs="i",yaxs="i")
points(zyg$position[zyg$contig=="3"], zyg$lnAOR[zyg$contig=="3"], pch=1, cex=0.5, col="green")
lines(xlim, rep( Aline,2), lty=3, col="red")
lines(xlim, rep(-Aline,2), lty=3, col="red")

plot(zyg$position[zyg$contig=="4"], zyg$lnBOR[zyg$contig=="4"], ylab="Chr4", pch=1, cex=0.5, col="red", xlim=xlim, ylim=ylim,  xaxt='n', xaxs="i",yaxs="i")
points(zyg$position[zyg$contig=="4"], zyg$lnAOR[zyg$contig=="4"], pch=1, cex=0.5, col="green")
lines(xlim, rep( Bline,2), lty=3, col="red")
lines(xlim, rep(-Bline,2), lty=3, col="red")

plot(zyg$position[zyg$contig=="5"], zyg$lnBOR[zyg$contig=="5"], ylab="Chr5", pch=1, cex=0.5, col="red", xlim=xlim, ylim=ylim,  xaxt='n', xaxs="i",yaxs="i")
points(zyg$position[zyg$contig=="5"], zyg$lnAOR[zyg$contig=="5"], pch=1, cex=0.5, col="green")
lines(xlim, rep( Bline,2), lty=3, col="red")
lines(xlim, rep(-Bline,2), lty=3, col="red")

plot(zyg$position[zyg$contig=="6"], zyg$lnBOR[zyg$contig=="6"], ylab="Chr6", pch=1, cex=0.5, col="red", xlim=xlim, ylim=ylim,  xaxt='n', xaxs="i",yaxs="i")
points(zyg$position[zyg$contig=="6"], zyg$lnAOR[zyg$contig=="6"], pch=1, cex=0.5, col="green")
lines(xlim, rep( Bline,2), lty=3, col="red")
lines(xlim, rep(-Bline,2), lty=3, col="red")

plot(zyg$position[zyg$contig=="7"], zyg$lnBOR[zyg$contig=="7"], ylab="Chr7", pch=1, cex=0.5, col="red", xlim=xlim, ylim=ylim,  xaxt='n', xaxs="i",yaxs="i")
points(zyg$position[zyg$contig=="7"], zyg$lnAOR[zyg$contig=="7"], pch=1, cex=0.5, col="green")
lines(xlim, rep( Bline,2), lty=3, col="red")
lines(xlim, rep(-Bline,2), lty=3, col="red")

plot(zyg$position[zyg$contig=="8"], zyg$lnBOR[zyg$contig=="8"], ylab="Chr8", pch=1, cex=0.5, col="red", xlim=xlim, ylim=ylim,  xaxt='n', xaxs="i",yaxs="i")
points(zyg$position[zyg$contig=="8"], zyg$lnAOR[zyg$contig=="8"], pch=1, cex=0.5, col="green")
lines(xlim, rep( Bline,2), lty=3, col="red")
lines(xlim, rep(-Bline,2), lty=3, col="red")

plot(zyg$position[zyg$contig=="9"], zyg$lnBOR[zyg$contig=="9"], ylab="Chr9", pch=1, cex=0.5, col="red", xlim=xlim, ylim=ylim,  xaxt='n', xaxs="i",yaxs="i")
points(zyg$position[zyg$contig=="9"], zyg$lnAOR[zyg$contig=="9"], pch=1, cex=0.5, col="green")
lines(xlim, rep( Bline,2), lty=3, col="red")
lines(xlim, rep(-Bline,2), lty=3, col="red")

plot(zyg$position[zyg$contig=="10"], zyg$lnBOR[zyg$contig=="10"], ylab="Chr10", pch=1, cex=0.5, col="red", xlim=xlim, ylim=ylim,  xaxt='n', xaxs="i",yaxs="i")
points(zyg$position[zyg$contig=="10"], zyg$lnAOR[zyg$contig=="10"], pch=1, cex=0.5, col="green")
lines(xlim, rep( Bline,2), lty=3, col="red")
lines(xlim, rep(-Bline,2), lty=3, col="red")

par(opar)
