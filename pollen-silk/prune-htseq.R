## prune out the genes that do not have W22/B73 snps from an htseq-counts file
## and write it back out

reference = "B73"
refSampleName = "SS354+355" # B73 silk REF
altSampleName = "PS422"     # W22 pollen ALT
mixSampleName = "S35464"    # B73 silk <- W22 pollen
mix1 = "S35464_1"
mix2 = "S35464_2"
mix3 = "S35464_4"
alt1 = "PS422_1"
alt2 = "PS422_2"
alt3 = "PS422_3"
ref1 = "SS354+355_0"
ref2 = "SS354+355_2"
ref3 = "SS354+355_3"

## reference = "W22"
## refSampleName = "SS364"     # W22 silk REF
## altSampleName = "YX24"      # B73 pollen ALT
## mixSampleName = "S364_354"  # W22 silk <- B73 pollen
## mix1 = "S364_354_0"
## mix2 = "S364_354_1"
## mix3 = "S364_354_3"
## alt1 = "YX24_2"
## alt2 = "YX24_5"
## alt3 = "YX24_10"
## ref1 = "SS364_1"
## ref2 = "SS364_2"
## ref3 = "SS364_3"

## reference = "B73"
## refSampleName = "YX24"      # B73 pollen REF
## altSampleName = "SS364"     # W22 silk ALT
## mixSampleName = "S364_354"  # B73 pollen -> W22 silk
## mix1 = "S364_354_0"
## mix2 = "S364_354_1"
## mix3 = "S364_354_3"
## alt1 = "SS364_1"
## alt2 = "SS364_2"
## alt3 = "SS364_3"
## ref1 = "YX24_10"
## ref2 = "YX24_2"
## ref3 = "YX24_5"

## reference = "W22"
## refSampleName = "PS422"     # W22 pollen REF
## altSampleName = "SS354+355" # B73 silk ALT
## mixSampleName = "S35464"    # W22 pollen -> B73 silk
## mix1 = "S35464_1"
## mix2 = "S35464_2"
## mix3 = "S35464_4"
## alt1 = "SS354+355_0"
## alt2 = "SS354+355_2"
## alt3 = "SS354+355_3"
## ref1 = "PS422_1"
## ref2 = "PS422_2"
## ref3 = "PS422_3"

## load the merged VCF records
merged = read.table(file=paste(reference,mixSampleName,refSampleName,"merged","txt",sep="."))

## we use the "clean" filter to remove marginal SNP calls
source("filter-merged.R")

## ## REF
## mix1.htseq = read.table(file=paste("STAR-",reference,"-",mix1,"/htseq-count.perfect.txt",sep=""), row.names=1, col.names=c("gene","count"))
## mix2.htseq = read.table(file=paste("STAR-",reference,"-",mix2,"/htseq-count.perfect.txt",sep=""), row.names=1, col.names=c("gene","count"))
## mix3.htseq = read.table(file=paste("STAR-",reference,"-",mix3,"/htseq-count.perfect.txt",sep=""), row.names=1, col.names=c("gene","count"))
## ref1.htseq = read.table(file=paste("STAR-",reference,"-",ref1,"/htseq-count.perfect.txt",sep=""), row.names=1, col.names=c("gene","count"))
## ref2.htseq = read.table(file=paste("STAR-",reference,"-",ref2,"/htseq-count.perfect.txt",sep=""), row.names=1, col.names=c("gene","count"))
## ref3.htseq = read.table(file=paste("STAR-",reference,"-",ref3,"/htseq-count.perfect.txt",sep=""), row.names=1, col.names=c("gene","count"))
## mix1.htseq.pruned = data.frame()
## mix2.htseq.pruned = data.frame()
## mix3.htseq.pruned = data.frame()
## ref1.htseq.pruned = data.frame()
## ref2.htseq.pruned = data.frame()
## ref3.htseq.pruned = data.frame()

## ALT
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
        ## BOTH
        mix1.htseq.pruned[gene,"count"] = mix1.htseq[gene,"count"]
        mix2.htseq.pruned[gene,"count"] = mix2.htseq[gene,"count"]
        mix3.htseq.pruned[gene,"count"] = mix3.htseq[gene,"count"]
        ## ## REF
        ## ref1.htseq.pruned[gene,"count"] = ref1.htseq[gene,"count"]
        ## ref2.htseq.pruned[gene,"count"] = ref2.htseq[gene,"count"]
        ## ref3.htseq.pruned[gene,"count"] = ref3.htseq[gene,"count"]
        ## ALT
        alt1.htseq.pruned[gene,"count"] = alt1.htseq[gene,"count"]
        alt2.htseq.pruned[gene,"count"] = alt2.htseq[gene,"count"]
        alt3.htseq.pruned[gene,"count"] = alt3.htseq[gene,"count"]
    } else {
        skipped = c(skipped, gene)
    }
}

## ## REF
## write.table(file=paste("STAR-",reference,"-",mix1,"/htseq-count.perfect.pruned.txt",sep=""), mix1.htseq.pruned, sep="\t", quote=F, col.names=F)
## write.table(file=paste("STAR-",reference,"-",mix2,"/htseq-count.perfect.pruned.txt",sep=""), mix2.htseq.pruned, sep="\t", quote=F, col.names=F)
## write.table(file=paste("STAR-",reference,"-",mix3,"/htseq-count.perfect.pruned.txt",sep=""), mix3.htseq.pruned, sep="\t", quote=F, col.names=F)
## write.table(file=paste("STAR-",reference,"-",ref1,"/htseq-count.perfect.pruned.txt",sep=""), ref1.htseq.pruned, sep="\t", quote=F, col.names=F)
## write.table(file=paste("STAR-",reference,"-",ref2,"/htseq-count.perfect.pruned.txt",sep=""), ref2.htseq.pruned, sep="\t", quote=F, col.names=F)
## write.table(file=paste("STAR-",reference,"-",ref3,"/htseq-count.perfect.pruned.txt",sep=""), ref3.htseq.pruned, sep="\t", quote=F, col.names=F)

## ALT
write.table(file=paste("STAR-",reference,"-",mix1,"/htseq-count.imperfect.pruned.txt",sep=""), mix1.htseq.pruned, sep="\t", quote=F, col.names=F)
write.table(file=paste("STAR-",reference,"-",mix2,"/htseq-count.imperfect.pruned.txt",sep=""), mix2.htseq.pruned, sep="\t", quote=F, col.names=F)
write.table(file=paste("STAR-",reference,"-",mix3,"/htseq-count.imperfect.pruned.txt",sep=""), mix3.htseq.pruned, sep="\t", quote=F, col.names=F)
write.table(file=paste("STAR-",reference,"-",alt1,"/htseq-count.imperfect.pruned.txt",sep=""), alt1.htseq.pruned, sep="\t", quote=F, col.names=F)
write.table(file=paste("STAR-",reference,"-",alt2,"/htseq-count.imperfect.pruned.txt",sep=""), alt2.htseq.pruned, sep="\t", quote=F, col.names=F)
write.table(file=paste("STAR-",reference,"-",alt3,"/htseq-count.imperfect.pruned.txt",sep=""), alt3.htseq.pruned, sep="\t", quote=F, col.names=F)

