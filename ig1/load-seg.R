##
## load the segregation data from a seg file
##                              ---------- counts ----------    Fisher's two-sided p
## contig start	REF	ALT	ref1	alt1	ref2	alt2	p
## 1	  31205	C	T	62	40	56	62	0.057868739309661944
## 1	  31213	T	C	44	55	39	68	0.25823090544022403
## 1	  31274	A	G	40	37	38	41	0.748876720004845
## 1	  43492	T	G	162	94	108	102	0.010963768277262271
##

segFile = readline(prompt="seg file: ")
seg = read.table(file=segFile, header=TRUE)
seg = seg[!is.nan(seg$p),]

## -log10(p)
seg$mlog10p = -log10(seg$p)

## odds ratio
seg$log10OR = log10((as.double(seg$ref1)*as.double(seg$alt2))/(as.double(seg$alt1)*as.double(seg$ref2)))

## log2(ALT/REF) -- useful for discerning true HETs
seg$logRatio1 = log2(seg$alt1/seg$ref1)
seg$logRatio2 = log2(seg$alt2/seg$ref2)

## ## comment out for chromosome-based seg file
## seg$tstart = seg$start
## seg$transcript = seg$contig
## for (i in 1:nrow(seg)) {
##     seqid = transcripts[seg$transcript[i],"seqid"]
##     start = transcripts[seg$transcript[i],"start"]
##     if (is.na(seqid)) {
##         seg$contig[i] = NA
##         seg$start[i] = NA
##     } else {
##         seg$contig[i] = seqid
##         seg$start[i] = start + seg$tstart[i] - 1
##     }
## }
## seg = seg[!is.na(seg$start),]


