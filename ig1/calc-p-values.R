##
## Add Fisher exact test p values to seg dataframe for Het/Hom, Het/NC, Hom/NC vs Non-Twin/Twin segregation
##

## Fisher's exact test p-value for the contingency table: Het/Hom Non-Twin/Twin OR Het/Ref OR Hom/Ref
for (i in 1:length(seg$pos)) {
    ## Het/Hom
    seg$pHetHom[i] = as.numeric(fisher.test(matrix(c(seg$Ahet[i],seg$Ahom[i],seg$Bhet[i],seg$Bhom[i]), nrow=2))["p.value"])
    ## Het/NC
    seg$pHetNC[i] = as.numeric(fisher.test(matrix(c(seg$Ahet[i],seg$Anc[i],seg$Bhet[i],seg$Bnc[i]), nrow=2))["p.value"])
    ## Hom/NC
    seg$pHomNC[i] = as.numeric(fisher.test(matrix(c(seg$Ahom[i],seg$Anc[i],seg$Bhom[i],seg$Bnc[i]), nrow=2))["p.value"])
    ## output if significant
    if (seg$pHetHom[i]<1e-2 || seg$pHetNC[i]<1e-2 || seg$pHomNC[i]<1e-2) {
        print(paste(seg$contig[i],seg$pos[i],
                    seg$Anc[i],seg$Ahet[i],seg$Ahom[i],
                    seg$Bnc[i],seg$Bhet[i],seg$Bhom[i],
                    seg$pHetNC[i],seg$pHomNC[i],seg$pHetHom[i]),
              quote=F)
    }
}
