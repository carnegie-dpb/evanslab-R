##
## plot -log10(p) calculated by SnpFilt caseControl in stacked plots
##

opar = par(mfrow=c(10,1))
par(mar=c(0.4,4,0.4,0.4))

maxLen = 305861025

xlim = c(0,maxLen)
yline1 = rep(-log10(0.05), 2)
yline2 = rep(-log10(0.01), 2)

plot(cc$Pos[cc$Chr=="1"], -log10(cc$CC_DOM[cc$Chr=="1"]), ylab="Chr1", pch=1, xlim=xlim,  xaxt='n')
lines(xlim, yline1, col="blue", lt=2)
lines(xlim, yline2, col="red", lt=2)

plot(cc$Pos[cc$Chr=="2"], -log10(cc$CC_DOM[cc$Chr=="2"]), ylab="Chr2", pch=1, xlim=xlim,  xaxt='n')
lines(xlim, yline1, col="blue", lt=2)
lines(xlim, yline2, col="red", lt=2)

plot(cc$Pos[cc$Chr=="3"], -log10(cc$CC_DOM[cc$Chr=="3"]), ylab="Chr3", pch=1, xlim=xlim,  xaxt='n')
lines(xlim, yline1, col="blue", lt=2)
lines(xlim, yline2, col="red", lt=2)

plot(cc$Pos[cc$Chr=="4"], -log10(cc$CC_DOM[cc$Chr=="4"]), ylab="Chr4", pch=1, xlim=xlim,  xaxt='n')
lines(xlim, yline1, col="blue", lt=2)
lines(xlim, yline2, col="red", lt=2)

plot(cc$Pos[cc$Chr=="5"], -log10(cc$CC_DOM[cc$Chr=="5"]), ylab="Chr5", pch=1, xlim=xlim,  xaxt='n')
lines(xlim, yline1, col="blue", lt=2)
lines(xlim, yline2, col="red", lt=2)

plot(cc$Pos[cc$Chr=="6"], -log10(cc$CC_DOM[cc$Chr=="6"]), ylab="Chr6", pch=1, xlim=xlim,  xaxt='n')
lines(xlim, yline1, col="blue", lt=2)
lines(xlim, yline2, col="red", lt=2)

plot(cc$Pos[cc$Chr=="7"], -log10(cc$CC_DOM[cc$Chr=="7"]), ylab="Chr7", pch=1, xlim=xlim,  xaxt='n')
lines(xlim, yline1, col="blue", lt=2)
lines(xlim, yline2, col="red", lt=2)

plot(cc$Pos[cc$Chr=="8"], -log10(cc$CC_DOM[cc$Chr=="8"]), ylab="Chr8", pch=1, xlim=xlim,  xaxt='n')
lines(xlim, yline1, col="blue", lt=2)
lines(xlim, yline2, col="red", lt=2)

plot(cc$Pos[cc$Chr=="9"], -log10(cc$CC_DOM[cc$Chr=="9"]), ylab="Chr9", pch=1, xlim=xlim,  xaxt='n')
lines(xlim, yline1, col="blue", lt=2)
lines(xlim, yline2, col="red", lt=2)

plot(cc$Pos[cc$Chr=="10"], -log10(cc$CC_DOM[cc$Chr=="10"]), ylab="Chr10", pch=1, xlim=xlim,  xaxt='n')
lines(xlim, yline1, col="blue", lt=2)
lines(xlim, yline2, col="red", lt=2)

par(opar)
