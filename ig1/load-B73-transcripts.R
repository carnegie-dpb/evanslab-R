library("ape")
library("taRifx")

transcripts = read.gff(file="~/genomes/Zm-B73-REFERENCE-GRAMENE-4.0/Zm00001d.2.gff3", GFF3=TRUE)

## get rid of factors
transcripts = remove.factors(transcripts)

## transcripts only!
transcripts = transcripts[transcripts$type=="mRNA",]
transcripts$length = transcripts$end - transcripts$start

## put the transcript IDs as row names
for (i in 1:length(transcripts$attributes)) {
    transcripts$ID[i] = substring(strsplit(transcripts$attributes[i], ";")[[1]][1], 15)
}
rownames(transcripts) = transcripts$ID
transcripts$ID = NULL

