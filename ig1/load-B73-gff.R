library("ape")

gff = read.gff(file="~/genomes/Zm-B73-REFERENCE-GRAMENE-4.0/Zm00001d.2.gff3", GFF3=TRUE)

## genes only!
gff = gff[gff$type=="gene",]

gff$length = gff$end - gff$start

## put the gene IDs as row names
for (i in 1:length(gff$attributes)) {
    gff$ID[i] = substring(strsplit(gff$attributes[i], ";")[[1]][1], 9)
}
rownames(gff) = gff$ID
gff$ID = NULL

