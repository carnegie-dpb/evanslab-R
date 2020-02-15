##
## load gene data from a genes-only GFF3 file
##

library(ape)

gffFile = "/home/shokin/genomes/Zm-B73-REFERENCE-GRAMENE-4.0/Zm00001d.2.gff3"

## genes only!
genes = read.gff(gffFile, GFF3=TRUE)
genes = genes[genes$type=="gene",]
genes$length = genes$end - genes$start

## put the gene IDs as row names
for (i in 1:length(genes$attributes)) {
    genes$ID[i] = substring(strsplit(genes$attributes[i], ";")[[1]][1], 9)
}
rownames(genes) = genes$ID
genes$ID = NULL
