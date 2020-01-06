##
## plot log Fisher p value in given region
##

plot.p = function(chr, start, end, geneIds=FALSE, boxcar=0) {
    
    strain = "B73"

    minsize = 100
    redline = 1e-6

    ylim = c(0,10)
    yline1 = rep(-log10(1e-2), 2)
    yline2 = rep(-log10(redline), 2)

    plotrange = seg$contig==chr & seg$start>=start & seg$start<=end & seg$size>=minsize
    highrange = seg$contig==chr & seg$start>=start & seg$start<=end & seg$size>=minsize & seg$p<redline

    plotseg = seg[plotrange,]

    if (boxcar>0) {
        ## add up the counts and THEN take Fisher's exact test
        mlog10p = rep(-1,boxcar) # don't show
        imin = boxcar + 1
        imax = nrow(plotseg) - boxcar
        for (i in imin:imax) {
            ilims = (i-boxcar):(i+boxcar)
            a = sum(plotseg$a[ilims])
            b = sum(plotseg$b[ilims])
            c = sum(plotseg$c[ilims])
            d = sum(plotseg$d[ilims])
            if ((a+b)>minsize && (c+d)>minsize) {
                p = fisher.test(x=matrix(c(a,b,c,d),nrow=2,ncol=2))$p.value
                mlog10p = c(mlog10p, -log10(p))
            } else {
                mlog10p = -1
            }
        }
        mlog10p = c(mlog10p, rep(-1,boxcar)) # don't show
    } else {
        mlog10p = plotseg$mlog10p
    }
    
    ## one plot for the given chromosome and region
    if (boxcar>0) {
        lines(plotseg$start, mlog10p, col="black", ylim=ylim,
             xlab=paste("Chr",chr,"Position"), ylab="-log10(p)")
    } else {
        plot(plotseg$start, mlog10p, pch=1, cex=0.5, col="black", ylim=ylim,
             xlab=paste("Chr",chr,"Position"), ylab="-log10(p)")
    }

    title(main=paste(strain,"calls: Fisher's Exact Test -log10(p) for ALT/REF counts vs Twin/Non-Twin\nCalls with at least",minsize,"total reads"), cex.main=0.9)
    
    if (!boxcar) {
        ## highlight highly significant p values with red
        points(seg$start[highrange], seg$mlog10p[highrange], pch=19, cex=0.5, col="darkred")
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
            lines(c(genes$start[i],genes$end[i]), c(0,0), col="blue", lwd=2)
            if (!boxcar) {
                generange = plotseg$start>=genes$start[i] & plotseg$start<=genes$end[i]
                points(plotseg$start[generange], plotseg$mlog10p[generange], pch=19, cex=0.1, col="blue")
            }
            if (geneIds) {
                geneId = sub("ID=gene:", "", strsplit(genes$attributes[i],";")[[1]][1])
                text(genes$start[i], -0.2, geneId, cex=0.8, pos=4, col="blue", offset=0)
            }
        }
    }
}

