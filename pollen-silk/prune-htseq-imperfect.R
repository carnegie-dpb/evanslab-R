## prune out the genes that do not have W22/B73 snps from an htseq-counts file
## and write it back out
##
## NOTE: you must have already loaded the merged VCF records and defined the experiment with read-merged.R

if (mixSampleName=="S35464" && refSampleName=="SS354+355" && altSampleName=="PS422") {
    mix1 = "S35464_1"
    mix2 = "S35464_2"
    mix3 = "S35464_4"
    alt1 = "PS422_1"
    alt2 = "PS422_2"
    alt3 = "PS422_3"
}

if (mixSampleName=="S364_354" && refSampleName=="SS364" && altSampleName=="YX24") {
    mix1 = "S364_354_0"
    mix2 = "S364_354_1"
    mix3 = "S364_354_3"
    alt1 = "YX24_2"
    alt2 = "YX24_5"
    alt3 = "YX24_10"
}

if (mixSampleName=="S364_354" && refSampleName=="YX24" && altSampleName=="SS364") {
    mix1 = "S364_354_0"
    mix2 = "S364_354_1"
    mix3 = "S364_354_3"
    alt1 = "SS364_1"
    alt2 = "SS364_2"
    alt3 = "SS364_3"
}

if (mixSampleName=="S35464" && refSampleName=="PS422" && altSampleName=="SS354+355") {
    mix1 = "S35464_1"
    mix2 = "S35464_2"
    mix3 = "S35464_4"
    alt1 = "SS354+355_0"
    alt2 = "SS354+355_2"
    alt3 = "SS354+355_3"
}

## we use the "clean" filter to remove marginal SNP calls
source("filter-merged.R")

mix1.htseq = read.table(file=paste("STAR-",reference,"-",mix1,"/htseq-count.imperfect.txt",sep=""), row.names=1, col.names=c("gene","count"))
mix2.htseq = read.table(file=paste("STAR-",reference,"-",mix2,"/htseq-count.imperfect.txt",sep=""), row.names=1, col.names=c("gene","count"))
mix3.htseq = read.table(file=paste("STAR-",reference,"-",mix3,"/htseq-count.imperfect.txt",sep=""), row.names=1, col.names=c("gene","count"))
alt1.htseq = read.table(file=paste("STAR-",reference,"-",alt1,"/htseq-count.imperfect.txt",sep=""), row.names=1, col.names=c("gene","count"))
alt2.htseq = read.table(file=paste("STAR-",reference,"-",alt2,"/htseq-count.imperfect.txt",sep=""), row.names=1, col.names=c("gene","count"))
alt3.htseq = read.table(file=paste("STAR-",reference,"-",alt3,"/htseq-count.imperfect.txt",sep=""), row.names=1, col.names=c("gene","count"))
mix1.htseq.pruned = data.frame()
mix2.htseq.pruned = data.frame()
mix3.htseq.pruned = data.frame()
alt1.htseq.pruned = data.frame()
alt2.htseq.pruned = data.frame()
alt3.htseq.pruned = data.frame()

skipped = c()

for (gene in rownames(mix1.htseq)) {
    if (length(merged$GENE[merged$GENE==gene & clean])>0) {
        mix1.htseq.pruned[gene,"count"] = mix1.htseq[gene,"count"]
        mix2.htseq.pruned[gene,"count"] = mix2.htseq[gene,"count"]
        mix3.htseq.pruned[gene,"count"] = mix3.htseq[gene,"count"]
        alt1.htseq.pruned[gene,"count"] = alt1.htseq[gene,"count"]
        alt2.htseq.pruned[gene,"count"] = alt2.htseq[gene,"count"]
        alt3.htseq.pruned[gene,"count"] = alt3.htseq[gene,"count"]
    } else {
        skipped = c(skipped, gene)
    }
}

write.table(file=paste("STAR-",reference,"-",mix1,"/htseq-count.imperfect.pruned.txt",sep=""), mix1.htseq.pruned, sep="\t", quote=F, col.names=F)
write.table(file=paste("STAR-",reference,"-",mix2,"/htseq-count.imperfect.pruned.txt",sep=""), mix2.htseq.pruned, sep="\t", quote=F, col.names=F)
write.table(file=paste("STAR-",reference,"-",mix3,"/htseq-count.imperfect.pruned.txt",sep=""), mix3.htseq.pruned, sep="\t", quote=F, col.names=F)
write.table(file=paste("STAR-",reference,"-",alt1,"/htseq-count.imperfect.pruned.txt",sep=""), alt1.htseq.pruned, sep="\t", quote=F, col.names=F)
write.table(file=paste("STAR-",reference,"-",alt2,"/htseq-count.imperfect.pruned.txt",sep=""), alt2.htseq.pruned, sep="\t", quote=F, col.names=F)
write.table(file=paste("STAR-",reference,"-",alt3,"/htseq-count.imperfect.pruned.txt",sep=""), alt3.htseq.pruned, sep="\t", quote=F, col.names=F)

