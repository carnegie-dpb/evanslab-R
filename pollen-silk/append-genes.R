## append the gene name, start, end and length to merged loci where appropriate
## REQUIRES load-genes.R

for (i in 1:nrow(merged)) {
    if (reference=="B73") {
        matches = genes.B73[genes.B73$seqid==merged$CHR[i] & genes.B73$start<=merged$POS[i] & genes.B73$end>=merged$POS[i],]
    } else if (reference=="W22") {
        matches = genes.W22[genes.W22$seqid==merged$CHR[i] & genes.W22$start<=merged$POS[i] & genes.W22$end>=merged$POS[i],]
    }
    if (dim(matches)[1]>0) {
        merged[i,"GENE"] = rownames(matches)[1]
        merged[i,"START"] = matches$start[1]
        merged[i,"END"] = matches$end[1]
        merged[i,"LENGTH"] = matches$length[1]
    }
}

## remove loci w/o gene association
merged = merged[!is.na(merged$GENE),]

