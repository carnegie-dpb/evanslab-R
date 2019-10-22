## easy routine to show the merged values in the given range

merged.range = function(chr, start, end) {
    return(merged[merged$CHR==chr & merged$POS>=start & merged$POS<=end,])
}
