##
## Add Fisher exact test p values to vcfcounts dataframe for Het/Hom, Het/NC, Hom/NC vs Non-Twin/Twin segregation
##

## Fisher's exact test p-value for the contingency table: Het/Hom Non-Twin/Twin OR Het/Ref OR Hom/Ref
for (i in 1:length(vcfcounts$pos)) {
    ## Het/Hom
    vcfcounts$pHetHom[i] = as.numeric(fisher.test(matrix(c(vcfcounts$Ahet[i],vcfcounts$Ahom[i],vcfcounts$Bhet[i],vcfcounts$Bhom[i]), nrow=2))["p.value"])
    ## Het/NC
    vcfcounts$pHetNC[i] = as.numeric(fisher.test(matrix(c(vcfcounts$Ahet[i],vcfcounts$Anc[i],vcfcounts$Bhet[i],vcfcounts$Bnc[i]), nrow=2))["p.value"])
    ## Hom/NC
    vcfcounts$pHomNC[i] = as.numeric(fisher.test(matrix(c(vcfcounts$Ahom[i],vcfcounts$Anc[i],vcfcounts$Bhom[i],vcfcounts$Bnc[i]), nrow=2))["p.value"])
    ## output if significant
    if (vcfcounts$pHetHom[i]<1e-2 || vcfcounts$pHetNC[i]<1e-2 || vcfcounts$pHomNC[i]<1e-2) {
        print(paste(vcfcounts$contig[i],vcfcounts$pos[i],
                    vcfcounts$Anc[i],vcfcounts$Ahet[i],vcfcounts$Ahom[i],
                    vcfcounts$Bnc[i],vcfcounts$Bhet[i],vcfcounts$Bhom[i],
                    vcfcounts$pHetNC[i],vcfcounts$pHomNC[i],vcfcounts$pHetHom[i]),
              quote=F)
    }
}
