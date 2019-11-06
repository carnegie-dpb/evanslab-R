## load the GFOLD diff results
## GeneSymbol	GeneName	GFOLD(0.01)	E-FDR	log2fdc	1stRPKM	2ndRPKM


##    STAR-B73-SS354+355/Aligned.sortedByCoord.f2.sorted.perfect 259835633       1
##    STAR-B73-S35464/Aligned.sortedByCoord.f2.sorted.perfect    533763006       1.93614
reference = "B73"
sample1 = "SS354+355"
sample2 = "S35464" 
norm1 = 1
norm2 = 259835633/533763006*1.93614

## ##   STAR-B73-YX24/Aligned.sortedByCoord.f2.sorted.perfect      93289322        1
## ##   STAR-B73-S364_354/Aligned.sortedByCoord.f2.sorted.perfect  102718908       24.5636
## reference = "B73"
## sample1 = "YX24"
## sample2 = "S364_354"
## norm1 = 1
## norm2 = 93289322/102718908*24.5636


## ##    STAR-W22-PS422/Aligned.sortedByCoord.f2.sorted.perfect     343363442       1
## ##    STAR-W22-S35464/Aligned.sortedByCoord.f2.sorted.perfect    432506945       17.3719
## reference = "W22"
## sample1 = "PS422"
## sample2 = "S35464" 
## norm1 = 1
## norm2 = 343363442/432506945*17.3719

## ##    STAR-W22-SS364/Aligned.sortedByCoord.f2.sorted.perfect     339665784       1.6779
## ##    STAR-W22-S364_354/Aligned.sortedByCoord.f2.sorted.perfect  205738646       1
## ## reference = "W22"
## sample1 = "SS364"
## sample2 = "S364_354"
## norm1 = 1.6779
## norm2 = 339665784/205738646*1

prefix = paste(sample1,".",sample2,sep="")
diffFile = paste("gfold-",reference,"/",prefix,".diff",sep="")

gfold = read.table(file=diffFile, sep="\t", header=FALSE, quote="")
colnames(gfold) = c("GeneSymbol","GeneName","GFOLD","E.FDR","log2fdc","RPKM1","RPKM2")
rownames(gfold) = gfold$GeneSymbol
gfold$GeneSymbol = NULL
gfold$E.FDR = NULL

## keep Zm0000* gene models
gfold = gfold[substr(rownames(gfold),1,6)=="Zm0000",]

## ## keep nonzero GFOLD values
## gfold = gfold[gfold$GFOLD!=0,]

## ## keep zero GFOLD values
## gfold = gfold[gfold$GFOLD==0,]
