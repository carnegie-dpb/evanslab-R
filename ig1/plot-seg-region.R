##
## plot log Fisher p value of each seg call on each chromosome, stacked plots
##
    
strain = "B73"

## chr = 1
## loc = 289612071
## extend = 2000
## gene = c(289610277,289611735)
## genename = "Zm00001d034333"

## chr = 1
## loc = 248631355
## extend = 3000

## chr = 5
## loc = 219939901
## extend = 10000
## gene = c(219939330,219947234)
## genename = "Zm00001d018392"

## chr = 9
## loc = 10359639
## extend = 30000

chr = 3
loc  =198728836
extend = 100

loc2 = 198728808
loc3 = 198728884

## chr = 3
## loc = 201850767
## extend = 1000

minsize = 100
redline = 1e-6
ymax = 20

ylim = c(0,ymax)
yline1 = rep(-log10(1e-2), 2)
yline2 = rep(-log10(redline), 2)

xmin = loc - extend
xmax = loc + extend
xlim = c(xmin,xmax)

lowrange  = seg$contig==chr & seg$size>=minsize & seg$start>=xmin & seg$start<=xmax
highrange = seg$contig==chr & seg$p<redline & seg$size>=minsize & seg$start>=xmin & seg$start<=xmax

## one plot for the given chromosome and region

plot(seg$start[lowrange], seg$mlog10p[lowrange],
     xlab=paste("Chr",chr,"Position"), ylab="-log10(p)", 
     pch=1, cex=0.5, col="black", xlim=xlim, ylim=ylim)

title(main=paste(strain,"calls: Fisher's Exact Test -log10(p) for ALT/REF counts vs Twin/Non-Twin\nCalls with at least",minsize,"total reads"), cex.main=0.9)
    
## highlight highly significant p values
points(seg$start[highrange], seg$mlog10p[highrange], pch=19, cex=0.5, col="darkred")
lines(xlim,  yline1, col="gray", lty=2)
lines(xlim,  yline2, col="gray", lty=2)

## vertical line on our center
lines(c(loc,loc), ylim, col="darkred", lty=3)
text(loc, ymax, loc, col="darkred")

lines(c(loc2,loc2), ylim, col="darkred", lty=3)
text(loc2, ymax, loc2, col="darkred")

lines(c(loc3,loc3), ylim, col="darkred", lty=3)
text(loc3, ymax, loc3, col="darkred")

## show genes from read.gff dataframe
## seqid  source type  start    end score strand phase
##     1 gramene gene  44289  49837    NA      +  <NA>

genes = gff[gff$seqid==chr & gff$start>=xmin & gff$start<=xmax,]
for (i in 1:length(genes$type)) {
    lines(c(genes$start[i],genes$end[i]), c(0,0), col="darkblue", lwd=2)
}

arrows(gene[1], 5, gene[2], 5, col="darkblue", lwd=5)
text(mean(gene), c(5,5), genename, col="darkblue", pos=3)

