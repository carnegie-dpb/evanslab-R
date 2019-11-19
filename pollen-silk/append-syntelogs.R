## Append the other line's syntelogs to the merged records
## REQUIRES that append-genes.R be run first!
for (i in 1:nrow(merged)) {
    if (reference=="B73") {
        merged$syntelogs.W22[i] = genes.B73$syntelogs.W22[rownames(genes.B73)==merged$GENE[i]]
    } else if (reference=="W22") {
        merged$syntelogs.B73[i] = genes.W22$syntelogs.B73[rownames(genes.W22)==merged$GENE[i]]
    }
}
