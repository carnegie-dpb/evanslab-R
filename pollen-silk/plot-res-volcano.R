## Volcano plot of DESeq2 results
##                        baseMean      log2FoldChange             lfcSE              stat             pvalue              padj      snps
##                       <numeric>           <numeric>         <numeric>         <numeric>          <numeric>         <numeric> <numeric>
## Zm00001d000001 236373.542774438  -0.303529307769506 0.209377951491212 -1.44967178066142  0.147150068617523 0.382775712497471         4

if (matches=="imperfect") {
    title = paste(reference,"imperfect match reads:",mixSampleName,"(mix) /",altSampleName,"(alt)")
} else {
    title = paste(reference,"perfect match reads:",mixSampleName,"(mix) /",refSampleName,"(ref)")
}

plot(res$log2FoldChange, -log10(res$pvalue),
     xlab="log2(FC)", ylab="-log10(p-value)",
     xlim=c(-10,10),
     main=title,
     cex=0.5
     )

## highlight significant
lowPadj = res$padj<5e-2
posLogFC = res$log2FoldChange >  1.0
negLogFC = res$log2FoldChange < -1.0
posSignificant = lowPadj & posLogFC
negSignificant = lowPadj & negLogFC

points(res$log2FoldChange[posSignificant], -log10(res$pvalue[posSignificant]), cex=0.3, pch=19, col="darkgreen")
text(res$log2FoldChange[posSignificant],   -log10(res$pvalue[posSignificant]), rownames(res)[posSignificant], cex=0.5, pos=4, offset=0.2, col="darkgreen")

points(res$log2FoldChange[negSignificant], -log10(res$pvalue[negSignificant]), cex=0.3, pch=19, col="darkred")
text(res$log2FoldChange[negSignificant],   -log10(res$pvalue[negSignificant]), rownames(res)[negSignificant], cex=0.5, pos=2, offset=0.2, col="darkred")
