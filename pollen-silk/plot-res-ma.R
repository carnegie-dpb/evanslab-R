## MA plot of DESeq2 results
##                        baseMean      log2FoldChange             lfcSE              stat             pvalue              padj
##                       <numeric>           <numeric>         <numeric>         <numeric>          <numeric>         <numeric>
## Zm00001d000001 236373.542774438  -0.303529307769506 0.209377951491212 -1.44967178066142  0.147150068617523 0.382775712497471

plot(res$baseMean, res$log2FoldChange, log="x", cex=0.5,
     xlab="mean base counts", ylab="log2(FC)",
     main=paste(reference,"perfect match reads:",mixSampleName,"(mix) /",refSampleName,"(ref)")
     )

## highlight significant
lowPadj = res$padj<5e-2
hasSNPs = res$snps>0
posLogFC = res$log2FoldChange >  1.0
negLogFC = res$log2FoldChange < -1.0
posSignificant = lowPadj & hasSNPs & posLogFC
negSignificant = lowPadj & hasSNPs & negLogFC

points(res$baseMean[posSignificant], res$log2FoldChange[posSignificant], cex=0.3, pch=19, col="darkgreen")
points(res$baseMean[negSignificant], res$log2FoldChange[negSignificant], cex=0.3, pch=19, col="darkred")

text(res$baseMean[posSignificant], res$log2FoldChange[posSignificant], rownames(res)[posSignificant], cex=0.5, pos=4, offset=0.2, col="darkgreen")
text(res$baseMean[negSignificant], res$log2FoldChange[negSignificant], rownames(res)[negSignificant], cex=0.5, pos=4, offset=0.2, col="darkred")
