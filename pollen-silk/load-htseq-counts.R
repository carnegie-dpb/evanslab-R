##
## load and merge the htseq-counts files
##
## NOTE: requires load-genes.R to load gene data
##

## reads.txt describes the experiment
reads.B73 = read.table("reads.B73.txt", sep="\t", row.names=1, col.names=c("id","tissue","samplename"))
reads.W22 = read.table("reads.W22.txt", sep="\t", row.names=1, col.names=c("id","tissue","samplename"))

filenames.B73 = paste(rownames(reads.B73), "/htseq-count.txt", sep="")
filenames.W22 = paste(rownames(reads.W22), "/htseq-count.txt", sep="")

htseq.B73 = read.table(filenames.B73[1], row.names=1, col.names=c("gene", reads.B73$samplename[1]))
for (i in 2:6) {
    counts = read.table(filenames.B73[i], row.names=1, col.names=c("gene",reads.B73$samplename[i]))
    htseq.B73[,reads.B73$samplename[i]] = counts[,reads.B73$samplename[i]]
}

htseq.W22 = read.table(filenames.W22[1], row.names=1, col.names=c("gene", reads.W22$samplename[1]))
for (i in 2:6) {
    counts = read.table(filenames.W22[i], row.names=1, col.names=c("gene",reads.W22$samplename[i]))
    htseq.W22[,reads.W22$samplename[i]] = counts[,reads.W22$samplename[i]]
}

## purge the genes that have less than 100 total counts
rowsums.B73 = rowSums(htseq.B73)
htseq.B73 = htseq.B73[rowsums.B73>=100,]

rowsums.W22 = rowSums(htseq.W22)
htseq.W22 = htseq.W22[rowsums.W22>=100,]

## yank the generic rows
htseq.B73 = subset(htseq.B73, rownames(htseq.B73)!="__no_feature")
htseq.B73 = subset(htseq.B73, rownames(htseq.B73)!="__ambiguous")
htseq.B73 = subset(htseq.B73, rownames(htseq.B73)!="__alignment_not_unique")

htseq.W22 = subset(htseq.W22, rownames(htseq.W22)!="__no_feature")
htseq.W22 = subset(htseq.W22, rownames(htseq.W22)!="__ambiguous")
htseq.W22 = subset(htseq.W22, rownames(htseq.W22)!="__alignment_not_unique")

## yank the ENSRNA rows
htseq.B73 = htseq.B73[substr(rownames(htseq.B73),1,6)!="ENSRNA",]
htseq.W22 = htseq.W22[substr(rownames(htseq.W22),1,6)!="ENSRNA",]

## determine the TPM values for each sample

## P.B73.B73 P.W22.B73 S.B73.P.W22.B73 S.W22.P.B73.B73 S.B73.B73
for (j in 1:length(rownames(htseq.B73))) {
    ID = rownames(htseq.B73)[j]
    htseq.B73$P.B73.B73.tpm[j] = htseq.B73$P.B73.B73[j]/genes.B73[ID,"length"]
    htseq.B73$P.W22.B73.tpm[j] = htseq.B73$P.W22.B73[j]/genes.B73[ID,"length"]
    htseq.B73$S.B73.P.W22.B73.tpm[j] = htseq.B73$S.B73.P.W22.B73[j]/genes.B73[ID,"length"]
    htseq.B73$S.W22.P.B73.B73.tpm[j] = htseq.B73$S.W22.P.B73.B73[j]/genes.B73[ID,"length"]
    htseq.B73$S.B73.B73.tpm[j] = htseq.B73$S.B73.B73[j]/genes.B73[ID,"length"]
}
## remove genes that aren't in the GFF
htseq.B73 = htseq.B73[!is.na(htseq.B73$P.B73.B73.tpm),]
## normalize to TPM
htseq.B73$P.B73.B73.tpm = htseq.B73$P.B73.B73.tpm/sum(htseq.B73$P.B73.B73.tpm)*1e6
htseq.B73$P.W22.B73.tpm = htseq.B73$P.W22.B73.tpm/sum(htseq.B73$P.W22.B73.tpm)*1e6
htseq.B73$S.B73.P.W22.B73.tpm = htseq.B73$S.B73.P.W22.B73.tpm/sum(htseq.B73$S.B73.P.W22.B73.tpm)*1e6
htseq.B73$S.W22.P.B73.B73.tpm = htseq.B73$S.W22.P.B73.B73.tpm/sum(htseq.B73$S.W22.P.B73.B73.tpm)*1e6
htseq.B73$S.B73.B73.tpm = htseq.B73$S.B73.B73.tpm/sum(htseq.B73$S.B73.B73.tpm)*1e6

## P.B73.W22 P.W22.W22 S.B73.P.W22.W22 S.W22.P.B73.W22 S.B73.W22      
for (j in 1:length(rownames(htseq.W22))) {
    ID = rownames(htseq.W22)[j]
    htseq.W22$P.B73.W22.tpm[j] = htseq.W22$P.B73.W22[j]/genes.W22[ID,"length"]
    htseq.W22$P.W22.W22.tpm[j] = htseq.W22$P.W22.W22[j]/genes.W22[ID,"length"]
    htseq.W22$S.B73.P.W22.W22.tpm[j] = htseq.W22$S.B73.P.W22.W22[j]/genes.W22[ID,"length"]
    htseq.W22$S.W22.P.B73.W22.tpm[j] = htseq.W22$S.W22.P.B73.W22[j]/genes.W22[ID,"length"]
    htseq.W22$S.B73.W22.tpm[j] = htseq.W22$S.B73.W22[j]/genes.W22[ID,"length"]
}
## remove genes that aren't in the GFF
htseq.W22 = htseq.W22[!is.na(htseq.W22$P.W22.W22.tpm),]
## normalize to TPM
htseq.W22$P.B73.W22.tpm = htseq.W22$P.B73.W22.tpm/sum(htseq.W22$P.B73.W22.tpm)*1e6
htseq.W22$P.W22.W22.tpm = htseq.W22$P.W22.W22.tpm/sum(htseq.W22$P.W22.W22.tpm)*1e6
htseq.W22$S.B73.P.W22.W22.tpm = htseq.W22$S.B73.P.W22.W22.tpm/sum(htseq.W22$S.B73.P.W22.W22.tpm)*1e6
htseq.W22$S.W22.P.B73.W22.tpm = htseq.W22$S.W22.P.B73.W22.tpm/sum(htseq.W22$S.W22.P.B73.W22.tpm)*1e6
htseq.W22$S.B73.W22.tpm = htseq.W22$S.B73.W22.tpm/sum(htseq.W22$S.B73.W22.tpm)*1e6
