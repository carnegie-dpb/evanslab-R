## append the number of SNPs from merged to res records, indicating that those genes can be distinguished between ref and alt components
for (gene in rownames(res)) {
    res[gene,"snps"] = length(merged$GENE[clean & merged$GENE==gene])
}
