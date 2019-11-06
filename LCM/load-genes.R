##
## load gene data from a genes-only GFF3 file
##

library(ape)

gffFile = "/home/shokin/genomes/Zm-B73-REFERENCE-GRAMENE-4.0/Zm00001d.2.gene.gff3"

genes = read.gff(gffFile, GFF3=TRUE)

## parse the gene ID out of the attributes
## ID=gene:Zm00001d027230;biotype=protein_coding;description=Mitochondrial transcription termination factor family protein;gene_id=Zm00001d027230;logic_name=maker_gene

## strsplit(genes$attributes,";")[[20]][1]
## [1] "ID=gene:Zm00001d027258"

genes$length = genes$end - genes$start

for (i in 1:length(genes$attributes)) {
    genes$ID[i] = substring(strsplit(genes$attributes[i], ";")[[1]][1], 9)
}
rownames(genes) = genes$ID
genes$ID = NULL
