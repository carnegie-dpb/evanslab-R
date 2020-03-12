library("ape")
library("taRifx")

## read and get rid of factors
gff = read.gff(file="~/genomes/Zm-B73-REFERENCE-GRAMENE-4.0/Zm00001d.2.gff3", GFF3=TRUE)
gff = remove.factors(gff)

## genes
genes = gff[gff$type=="gene",]
genes$length = genes$end - genes$start

## CDSes
cds = gff[gff$type=="CDS",]
cds$length = cds$end - cds$start

## put the gene IDs as row names
## ID=gene:Zm00001d027230;biotype=protein_coding;description=Zm00001d027230;gene_id=Zm00001d027230;logic_name=maker_gene
for (i in 1:nrow(genes)) {
    genes$ID[i] = substring(strsplit(genes$attributes[i], ";")[[1]][1], 9)
}
rownames(genes) = genes$ID
genes$ID = NULL

## add the CDS IDs and gene names
## ID=CDS:Zm00001d027230_P001;Parent=transcript:Zm00001d027230_T001;protein_id=Zm00001d027230_P001
for (i in 1:nrow(cds)) {
    cds$ID[i] = substring(strsplit(cds$attributes[i], ";")[[1]][1], 8)
    cds$gene[i] = substring(cds$ID[i], 1, 14)
}

