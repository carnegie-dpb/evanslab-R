##
## load B73 gene data from a GFF3 file
##

library(ape)

gffFile.B73 = "/home/shokin/genomes/Zm-B73-REFERENCE-GRAMENE-4.0/Zm00001d.2.gff3"

## load GFF files
print("Loading GFF...", quote=F)
genes.B73 = read.gff(gffFile.B73, GFF3=TRUE)

## extract genes
print("Extracting genes...", quote=F)
genes.B73 = genes.B73[genes.B73$type=="gene",]

## get lengths
genes.B73$length = genes.B73$end - genes.B73$start

## drop scaffolds
genes.B73 = genes.B73[substr(genes.B73$seqid,1,9)!="B73V4_ctg",]

## parse the gene ID out of the B73 attributes
## ID=gene:Zm00001d027230;biotype=protein_coding;description=Mitochondrial transcription termination factor family protein;gene_id=Zm00001d027230;logic_name=maker_gene
## strsplit(genes.B73$attributes,";")[[20]][1]
## [1] "ID=gene:Zm00001d027258"
for (i in 1:nrow(genes.B73)) {
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

## now load the syntelogs for each line's genes; there can be multiple rows per gene
## B73.W22.syntelogs
## B73v4_gene_model_id B73v4_chr B73v4_start B73v4_stop B73v4_strand W22_gene_model_id W22_chr W22_start W22_stop W22_strand   X.ID   eval   score
## Zm00001d027230      1         44351       47995      1            Zm00004b000001    1       12652     16389    1            100.00 1e-250 50
print("Loading syntelogs...", quote=F)
source("load-B73-W22-synmap.R")
for (i in 1:nrow(genes.B73)) {
    syntelogs = B73.W22.syntelogs$W22_gene_model_id[B73.W22.syntelogs$B73v4_gene_model_id==rownames(genes.B73)[i]]
    if (length(syntelogs)>0) {
        genes.B73$syntelogs.W22[i] = paste(syntelogs, collapse=" ")
    }
}
