## Append the other line's syntelogs to the res records
for (gene in rownames(res)) {
    if (reference=="B73") {
        res[gene,"syntelogs.W22"] = genes.B73$syntelogs.W22[rownames(genes.B73)==gene]
    } else if (reference=="W22") {
        res[gene,"syntelogs.B73"] = genes.W22$syntelogs.B73[rownames(genes.W22)==gene]
    }
}
