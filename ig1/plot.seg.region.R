##
## plot log Fisher p value in given region
##

plot.seg.region = function(chr, start, end, geneIds=FALSE, boxcar=0) {
    
    strain = "B73"

    minsize = 100
    redline = 1e-6
    if (boxcar>0) {
        ymax = 5
    } else {
        ymax = 20
    }

    ylim = c(0,ymax)
    yline1 = rep(-log10(1e-2), 2)
    yline2 = rep(-log10(redline), 2)

    xlim = c(start,end)

    plotrange = seg$contig==chr & seg$start>=start & seg$start<=end
    lowrange  = seg$contig==chr & seg$start>=start & seg$start<=end & seg$size>=minsize
    highrange = seg$contig==chr & seg$start>=start & seg$start<=end & seg$size>=minsize & seg$p<redline

    plotseg = seg[plotrange,]

    if (boxcar>0) {
        ## boxcar average counts -> p
        mlog10p = rep(-1,boxcar)
        imin = boxcar + 1
        imax = nrow(plotseg) - boxcar
        for (i in imin:imax) {
            ilims = (i-boxcar):(i+boxcar)
            mlog10p = c(mlog10p, mean(plotseg$mlog10p[ilims]))
        }
        mlog10p = c(mlog10p, rep(-1,boxcar))
    } else {
        mlog10p = plotseg$mlog10p
    }
    
    ## one plot for the given chromosome and region
    plot(plotseg$start, mlog10p,
         xlab=paste("Chr",chr,"Position"), ylab="-log10(p)", 
         pch=1, cex=0.5, col="black", xlim=xlim, ylim=ylim)

    title(main=paste(strain,"calls: Fisher's Exact Test -log10(p) for ALT/REF counts vs Twin/Non-Twin\nCalls with at least",minsize,"total reads"), cex.main=0.9)
    
    if (!boxcar) {
        ## highlight highly significant p values
        points(seg$start[highrange], seg$mlog10p[highrange], pch=19, cex=0.5, col="darkred")
        lines(xlim,  yline1, col="gray", lty=2)
        lines(xlim,  yline2, col="red", lty=2)
    }

    ## show genes from read.gff dataframe
    ## seqid  source type  start    end score strand phase
    ##     1 gramene gene  44289  49837    NA      +  <NA>


    ## seqid  source type start   end score strand phase                                                 attributes
    ##     1 gramene gene 44289 49837    NA      +  <NA>  ID=gene:Zm00001d027230;biotype=protein_coding;description=Zm00001d027230;gene_id=Zm00001d027230;logic_name=maker_gene

    ## [[1]]
    ## [1] "ID=gene:Zm00001d027230"     "biotype=protein_coding"    
    ## [3] "description=Zm00001d027230" "gene_id=Zm00001d027230"    
    ## [5] "logic_name=maker_gene"

    genes = gff[gff$type=="gene" & gff$seqid==chr & gff$start>=start & gff$start<=end,]
    if (nrow(genes)>0) {
        for (i in 1:nrow(genes)) {
            lines(c(genes$start[i],genes$end[i]), c(0,0), col="blue", lwd=1)
            if (geneIds) {
                geneId = sub("ID=gene:", "", strsplit(genes$attributes[i],";")[[1]][1])
                text(genes$start[i], 0, geneId, cex=0.8, pos=2)
            }
        }
    }
}

