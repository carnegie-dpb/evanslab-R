##
## load and merge the htseq-counts files
##
## htseq/gsnap-B73-ACMLCM1.htseq-count.txt
##
## NOTE: requires load-genes.R to load gene data

## reads.txt describes the experiment
reads = read.table("reads.txt", sep="\t", row.names=1, col.names=c("id","tissue","samplename"))

## read and merge em'
filename = paste("htseq/gsnap-B73-",rownames(reads[1,]),".htseq-count.txt",sep="")
htseq.counts = read.table(filename, row.names=1, col.names=c("gene",reads$samplename[1]))
for (i in 2:length(reads$tissue)) {
    filename = paste("htseq/gsnap-B73-",rownames(reads[i,]),".htseq-count.txt",sep="")
    samplename = reads$samplename[i]
    counts = read.table(filename, row.names=1, col.names=c("gene",samplename))
    htseq.counts[,samplename] = counts[,samplename]
}

## purge the genes that have less than 20 total counts
rowsums = rowSums(htseq.counts)
htseq.counts = htseq.counts[rowsums>=20,]

## yank the generic rows
htseq.counts = subset(htseq.counts, rownames(htseq.counts)!="__no_feature")
htseq.counts = subset(htseq.counts, rownames(htseq.counts)!="__ambiguous")
htseq.counts = subset(htseq.counts, rownames(htseq.counts)!="__alignment_not_unique")

## yank the ENSRNA rows
htseq.counts = htseq.counts[substr(rownames(htseq.counts),1,6)!="ENSRNA",]

## determine the TPM values for each sample
for (j in 1:length(rownames(htseq.counts))) {
    ID = rownames(htseq.counts)[j]
    htseq.counts$WES1.tpm[j] = htseq.counts$WES1[j]/genes[ID,"length"]
    htseq.counts$WES2.tpm[j] = htseq.counts$WES2[j]/genes[ID,"length"]
    htseq.counts$WES3.tpm[j] = htseq.counts$WES3[j]/genes[ID,"length"]
    htseq.counts$WES4.tpm[j] = htseq.counts$WES4[j]/genes[ID,"length"]
    htseq.counts$ACC1.tpm[j] = htseq.counts$ACC1[j]/genes[ID,"length"]
    htseq.counts$ACC2.tpm[j] = htseq.counts$ACC2[j]/genes[ID,"length"]
    htseq.counts$ACC3.tpm[j] = htseq.counts$ACC3[j]/genes[ID,"length"]
    htseq.counts$ACC4.tpm[j] = htseq.counts$ACC4[j]/genes[ID,"length"]
    htseq.counts$CC1.tpm[j] = htseq.counts$CC1[j]/genes[ID,"length"]
    htseq.counts$CC2.tpm[j] = htseq.counts$CC2[j]/genes[ID,"length"]
    htseq.counts$CC3.tpm[j] = htseq.counts$CC3[j]/genes[ID,"length"]
    htseq.counts$CC4.tpm[j] = htseq.counts$CC4[j]/genes[ID,"length"]
    htseq.counts$EA1.tpm[j] = htseq.counts$EA1[j]/genes[ID,"length"]
    htseq.counts$EA2.tpm[j] = htseq.counts$EA2[j]/genes[ID,"length"]
    htseq.counts$EA3.tpm[j] = htseq.counts$EA3[j]/genes[ID,"length"]
    htseq.counts$EA4.tpm[j] = htseq.counts$EA4[j]/genes[ID,"length"]
    htseq.counts$N1.tpm[j] = htseq.counts$N1[j]/genes[ID,"length"]
    htseq.counts$N2.tpm[j] = htseq.counts$N2[j]/genes[ID,"length"]
    htseq.counts$N3.tpm[j] = htseq.counts$N3[j]/genes[ID,"length"]
    htseq.counts$N4.tpm[j] = htseq.counts$N4[j]/genes[ID,"length"]
    htseq.counts$P1.tpm[j] = htseq.counts$P1[j]/genes[ID,"length"]
    htseq.counts$P2.tpm[j] = htseq.counts$P2[j]/genes[ID,"length"]
    htseq.counts$P3.tpm[j] = htseq.counts$P3[j]/genes[ID,"length"]
    htseq.counts$P4.tpm[j] = htseq.counts$P4[j]/genes[ID,"length"]
}

## remove genes that aren't in the GFF
htseq.counts = htseq.counts[!is.na(htseq.counts$WES1.tpm),]

## normalize to TPM
htseq.counts$WES1.tpm = htseq.counts$WES1.tpm/sum(htseq.counts$WES1.tpm)*1e6
htseq.counts$WES2.tpm = htseq.counts$WES2.tpm/sum(htseq.counts$WES2.tpm)*1e6
htseq.counts$WES3.tpm = htseq.counts$WES3.tpm/sum(htseq.counts$WES3.tpm)*1e6
htseq.counts$WES4.tpm = htseq.counts$WES4.tpm/sum(htseq.counts$WES4.tpm)*1e6
htseq.counts$ACC1.tpm = htseq.counts$ACC1.tpm/sum(htseq.counts$ACC1.tpm)*1e6
htseq.counts$ACC2.tpm = htseq.counts$ACC2.tpm/sum(htseq.counts$ACC2.tpm)*1e6
htseq.counts$ACC3.tpm = htseq.counts$ACC3.tpm/sum(htseq.counts$ACC3.tpm)*1e6
htseq.counts$ACC4.tpm = htseq.counts$ACC4.tpm/sum(htseq.counts$ACC4.tpm)*1e6
htseq.counts$CC1.tpm = htseq.counts$CC1.tpm/sum(htseq.counts$CC1.tpm)*1e6
htseq.counts$CC2.tpm = htseq.counts$CC2.tpm/sum(htseq.counts$CC2.tpm)*1e6
htseq.counts$CC3.tpm = htseq.counts$CC3.tpm/sum(htseq.counts$CC3.tpm)*1e6
htseq.counts$CC4.tpm = htseq.counts$CC4.tpm/sum(htseq.counts$CC4.tpm)*1e6
htseq.counts$EA1.tpm = htseq.counts$EA1.tpm/sum(htseq.counts$EA1.tpm)*1e6
htseq.counts$EA2.tpm = htseq.counts$EA2.tpm/sum(htseq.counts$EA2.tpm)*1e6
htseq.counts$EA3.tpm = htseq.counts$EA3.tpm/sum(htseq.counts$EA3.tpm)*1e6
htseq.counts$EA4.tpm = htseq.counts$EA4.tpm/sum(htseq.counts$EA4.tpm)*1e6
htseq.counts$N1.tpm = htseq.counts$N1.tpm/sum(htseq.counts$N1.tpm)*1e6
htseq.counts$N2.tpm = htseq.counts$N2.tpm/sum(htseq.counts$N2.tpm)*1e6
htseq.counts$N3.tpm = htseq.counts$N3.tpm/sum(htseq.counts$N3.tpm)*1e6
htseq.counts$N4.tpm = htseq.counts$N4.tpm/sum(htseq.counts$N4.tpm)*1e6
htseq.counts$P1.tpm = htseq.counts$P1.tpm/sum(htseq.counts$P1.tpm)*1e6
htseq.counts$P2.tpm = htseq.counts$P2.tpm/sum(htseq.counts$P2.tpm)*1e6
htseq.counts$P3.tpm = htseq.counts$P3.tpm/sum(htseq.counts$P3.tpm)*1e6
htseq.counts$P4.tpm = htseq.counts$P4.tpm/sum(htseq.counts$P4.tpm)*1e6

## average the TPMs over the bioreps for each tissue
htseq.counts$WES.tpm = (htseq.counts$WES1.tpm + htseq.counts$WES2.tpm + htseq.counts$WES3.tpm + htseq.counts$WES4.tpm) / 4
htseq.counts$ACC.tpm = (htseq.counts$ACC1.tpm + htseq.counts$ACC2.tpm + htseq.counts$ACC3.tpm + htseq.counts$ACC4.tpm) / 4
htseq.counts$CC.tpm  = (htseq.counts$CC1.tpm  + htseq.counts$CC2.tpm  + htseq.counts$CC3.tpm  + htseq.counts$CC4.tpm) / 4
htseq.counts$EA.tpm  = (htseq.counts$EA1.tpm  + htseq.counts$EA2.tpm  + htseq.counts$EA3.tpm  + htseq.counts$EA4.tpm) / 4
htseq.counts$N.tpm   = (htseq.counts$N1.tpm   + htseq.counts$N2.tpm   + htseq.counts$N3.tpm   + htseq.counts$N4.tpm) / 4
htseq.counts$P.tpm   = (htseq.counts$P1.tpm   + htseq.counts$P2.tpm   + htseq.counts$P3.tpm   + htseq.counts$P4.tpm) / 4

## create the tissue total percentage values for fractions
## NOTE: WES includes ACC, CC and EA contributions, so do NOT include it in the total for those fractions!
## WES can be compared to N and P.
htseq.counts$NonWES = htseq.counts$ACC.tpm + htseq.counts$CC.tpm + htseq.counts$EA.tpm + htseq.counts$N.tpm + htseq.counts$P.tpm
htseq.counts$WES.N.P = htseq.counts$WES.tpm + htseq.counts$N.tpm + htseq.counts$P.tpm

## calculate the appropriate fractions
htseq.counts$WES.frac = htseq.counts$WES.tpm / htseq.counts$WES.N.P
htseq.counts$ACC.frac = htseq.counts$ACC.tpm / htseq.counts$NonWES
htseq.counts$CC.frac  = htseq.counts$CC.tpm / htseq.counts$NonWES
htseq.counts$EA.frac  = htseq.counts$EA.tpm / htseq.counts$NonWES
htseq.counts$N.frac   = htseq.counts$N.tpm / htseq.counts$NonWES
htseq.counts$P.frac   = htseq.counts$P.tpm / htseq.counts$NonWES
