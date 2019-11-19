##
## load B73 gene data from a GFF3 file
##

library(ape)

gffFile.W22 = "/home/shokin/genomes/Zm-W22-REFERENCE-NRGENE-2.0/Zm00004b.gff"

## load GFF files
print("Loading GFF...", quote=F)
genes.W22 = read.gff(gffFile.W22, GFF3=TRUE)

## extract genes
print("Extracting genes...", quote=F)
genes.W22 = genes.W22[genes.W22$type=="gene",]

## get lengths
genes.W22$length = genes.W22$end - genes.W22$start

## drop scaffolds
genes.W22 = genes.W22[substr(genes.W22$seqid,1,8)!="scaffold",]

## parse the gene ID out of the W22 attributes
## ID=Zm00004b000001;Name=Zm00004b000001;biotype=protein_coding
## strsplit(genes.W22$attributes,";")[[20]][1]
## [1] "ID=Zm00004b000020"
for (i in 1:nrow(genes.W22)) {
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

## now load the syntelogs from the other line; there can be multiple rows per gene
## B73.W22.syntelogs
## B73v4_gene_model_id B73v4_chr B73v4_start B73v4_stop B73v4_strand W22_gene_model_id W22_chr W22_start W22_stop W22_strand   X.ID   eval   score
## Zm00001d027230      1         44351       47995      1            Zm00004b000001    1       12652     16389    1            100.00 1e-250 50
print("Loading syntelogs...", quote=F)
source("load-B73-W22-synmap.R")
for (i in 1:nrow(genes.W22)) {
    syntelogs = B73.W22.syntelogs$B73v4_gene_model_id[B73.W22.syntelogs$W22_gene_model_id==rownames(genes.W22)[i]]
    if (length(syntelogs)>0) {
        genes.W22$syntelogs.B73[i] = paste(syntelogs, collapse=" ")
    }
}
