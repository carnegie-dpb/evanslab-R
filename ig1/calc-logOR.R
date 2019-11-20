##
## add log odds-ratio columns to seg dataframe
##

## need totals for odds-ratio stats
A = 45
B = 48

## HET log odds ratio: restrict to cases where we have at least 10 calls from each group
for (i in 1:length(seg$Ahet)) {
    if (seg$Ahet[i]>=10 && seg$Bhet[i]>=10) {
        seg$logHetOR = log( seg$Bhet/(B-seg$Bhet) / ( seg$Ahet/(A-seg$Ahet) ) )
        print(paste(seg$contig[i],seg$pos[i],seg$logHetOR[i]))
    } else {
        seg$logHetOR = 0.0
    }
}


