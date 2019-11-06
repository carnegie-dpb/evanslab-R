##
## load gene data from a GFF3 file
##

library(ape)

gffFile.B73 = "/home/shokin/genomes/Zm-B73-REFERENCE-GRAMENE-4.0/Zm00001d.2.gff3"
gffFile.W22 = "/home/shokin/genomes/Zm-W22-REFERENCE-NRGENE-2.0/Zm00004b.gff"

genes.B73 = read.gff(gffFile.B73, GFF3=TRUE)
genes.W22 = read.gff(gffFile.W22, GFF3=TRUE)

## extract genes
genes.B73 = genes.B73[genes.B73$type=="gene",]
genes.W22 = genes.W22[genes.W22$type=="gene",]

## get lengths
genes.B73$length = genes.B73$end - genes.B73$start
genes.W22$length = genes.W22$end - genes.W22$start

## drop scaffolds
genes.W22 = genes.W22[substr(genes.W22$seqid,1,8)!="scaffold",]
genes.B73 = genes.B73[substr(genes.B73$seqid,1,9)!="B73V4_ctg",]

## parse the gene ID out of the B73 attributes
## ID=gene:Zm00001d027230;biotype=protein_coding;description=Mitochondrial transcription termination factor family protein;gene_id=Zm00001d027230;logic_name=maker_gene
## strsplit(genes.B73$attributes,";")[[20]][1]
## [1] "ID=gene:Zm00001d027258"
for (i in 1:length(genes.B73$attributes)) {
    genes.B73$ID[i] = substring(strsplit(genes.B73$attributes[i], ";")[[1]][1], 9)
}
rownames(genes.B73) = genes.B73$ID
genes.B73$ID = NULL
genes.B73$source = NULL
genes.B73$type = NULL
genes.B73$score = NULL
genes.B73$strand = NULL
genes.B73$phase = NULL
genes.B73$attributes = NULL

## parse the gene ID out of the W22 attributes
## ID=Zm00004b000001;Name=Zm00004b000001;biotype=protein_coding
## strsplit(genes.W22$attributes,";")[[20]][1]
## [1] "ID=Zm00004b000020"
for (i in 1:length(genes.W22$attributes)) {
    genes.W22$ID[i] = substring(strsplit(genes.W22$attributes[i], ";")[[1]][1], 4)
}
rownames(genes.W22) = genes.W22$ID
genes.W22$ID = NULL
genes.W22$source = NULL
genes.W22$type = NULL
genes.W22$score = NULL
genes.W22$strand = NULL
genes.W22$phase = NULL
genes.W22$attributes = NULL
