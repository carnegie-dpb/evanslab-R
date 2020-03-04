##
## load the segregation data from a seg file
##
## a=NON-TWIN REF
## b=NON-TWIN ALT
## c=TWIN REF
## d=TWIN ALT
##
## NOTE: there may be consecutive calls at the same position, e.g.
## contig                  start                                       REF ALT   a  b  c  d size          p     mlog10p signif
##                   1 102591307                                         T   C  11 50  5 34  122 0.17886000   0.7474868  false
##                   1 102591307 TAGCCAGGGCTTGGTCTCGGGTGAGTCGTGACTATGTCTCC   T  93 44 75 29  274 0.08815873   1.0547347  false
##   GRMZM5G800096_T01       295 C                                           G   1 63  0 76  140 0.45402884   0.3429166  false
## Zm00001d000037_T001       997 G                                           C 114 26 95 52 280  6.212756e-04 3.2067157   true

segFile = readline(prompt="seg file: ")
seg = read.table(file=segFile, header=TRUE)
seg = seg[!is.nan(seg$p),]
seg$log10OR = log10((as.double(seg$a)*as.double(seg$d))/(as.double(seg$b)*as.double(seg$c)))

## comment out for chromosome-based seg file
seg$tstart = seg$start
seg$transcript = seg$contig
for (i in 1:nrow(seg)) {
    seqid = transcripts[seg$transcript[i],"seqid"]
    start = transcripts[seg$transcript[i],"start"]
    if (is.na(seqid)) {
        seg$contig[i] = NA
        seg$start[i] = NA
    } else {
        seg$contig[i] = seqid
        seg$start[i] = start + seg$tstart[i] - 1
    }
}

## cull the non-identified transcripts
seg = seg[!is.na(seg$start),]


