##
## purge outlier genes from htseq.counts
##

htseq.counts = subset(htseq.counts, rownames(htseq.counts)!="Zm00001d000189")
htseq.counts = subset(htseq.counts, rownames(htseq.counts)!="Zm00001d000477")
htseq.counts = subset(htseq.counts, rownames(htseq.counts)!="Zm00001d005080")
htseq.counts = subset(htseq.counts, rownames(htseq.counts)!="Zm00001d011014")
htseq.counts = subset(htseq.counts, rownames(htseq.counts)!="Zm00001d011017")
htseq.counts = subset(htseq.counts, rownames(htseq.counts)!="Zm00001d011026")
htseq.counts = subset(htseq.counts, rownames(htseq.counts)!="Zm00001d011028")
htseq.counts = subset(htseq.counts, rownames(htseq.counts)!="Zm00001d018252")
htseq.counts = subset(htseq.counts, rownames(htseq.counts)!="Zm00001d019605")
htseq.counts = subset(htseq.counts, rownames(htseq.counts)!="Zm00001d023090")
htseq.counts = subset(htseq.counts, rownames(htseq.counts)!="Zm00001d029056")
htseq.counts = subset(htseq.counts, rownames(htseq.counts)!="Zm00001d037319")
htseq.counts = subset(htseq.counts, rownames(htseq.counts)!="Zm00001d046613")

htseq.counts = subset(htseq.counts, rownames(htseq.counts)!="GRMZM5G806488")
htseq.counts = subset(htseq.counts, rownames(htseq.counts)!="GRMZM5G807592")
htseq.counts = subset(htseq.counts, rownames(htseq.counts)!="GRMZM5G844030")
