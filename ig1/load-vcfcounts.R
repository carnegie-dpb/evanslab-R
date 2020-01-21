## load the vcfcounts GBS data
## chr     pos     ref     alt     twHet   twHom   twNC    ntHet   ntHom   ntNC
## 1       33879   C*      T       1       3       44      0       3       42


vcfcounts = read.table(file="vcfcounts.txt", header=TRUE)

vcfcounts$log10OR = log10((vcfcounts$twHet/48)/(vcfcounts$ntHet/45))

## for (i in 1:nrow(vcfcounts)) {
##     contingency = matrix(c(vcfcounts$ntHet[i],vcfcounts$twHet[i],45,48),nrow=2,ncol=2)
##     vcfcounts$pHet[i] = fisher.test(contingency)$p.value
##     if (vcfcounts$pHet[i]<1e-2) {
##         print(vcfcounts[i,])
##     }
## }
