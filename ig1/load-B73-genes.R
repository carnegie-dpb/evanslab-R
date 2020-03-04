library("ape")
library("taRifx")

genes = read.gff(file="~/genomes/Zm-B73-REFERENCE-GRAMENE-4.0/Zm00001d.2.gff3", GFF3=TRUE)

## get rid of factors
genes = remove.factors(genes)

## genes only!
genes = genes[genes$type=="gene",]

genes$length = genes$end - genes$start

## put the gene IDs as row names
for (i in 1:length(genes$attributes)) {
    genes$ID[i] = substring(strsplit(genes$attributes[i], ";")[[1]][1], 9)
}
rownames(genes) = genes$ID
genes$ID = NULL

