source("plot.logp.R")

##
## plot a series of PNG files with increasing values of minHetSamples for the given value of minReads
##
seg.movie = function(minReads) {

    for (i in 1:18*5) {
        frame = i
        if (i<10) frame = paste("0",i,sep="")
        fileName = paste("plot.seg.",minReads,".",frame,".png", sep="")
        print(fileName)
        png(fileName, width=1600, height=900)
        plot.logp(minReads, minHetSamples=i)
        dev.off()
    }

}
    
