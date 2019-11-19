## merge the data frames that are loaded with loadVCF
source("loadVCF.R")

singleTissues = c(
    "SS354+355"="B73 silk", 
    "YX24"="B73 pollen", 
    "SS364"="W22 silk", 
    "PS422"="W22 pollen"
)

## reference = "B73"
## refMixRatio = 1.81769
## altMixRatio = 0.04887
## mixSampleName = "S35464"    # B73 silk <- W22 pollen
## refSampleName = "SS354+355" # B73 silk REF
## altSampleName = "PS422"     # W22 pollen ALT

## reference = "B73"
## refMixRatio = 0.23941
## altMixRatio = 0.55649
## mixSampleName = "S364_354"  # B73 pollen -> W22 silk
## refSampleName = "YX24"      # B73 pollen REF
## altSampleName = "SS364"     # W22 silk ALT

## reference = "W22"
## refMixRatio = 0.57338
## altMixRatio = 0.22939
## mixSampleName = "S364_354"  # W22 silk <- B73 pollen
## refSampleName = "SS364"     # W22 silk REF
## altSampleName = "YX24"      # B73 pollen ALT

## reference = "W22"
## refMixRatio = 0.04977
## altMixRatio = 1.81575
## mixSampleName = "S35464"    # W22 pollen -> B73 silk
## refSampleName = "PS422"     # W22 pollen REF
## altSampleName = "SS354+355" # B73 silk ALT

print(paste("Reference=",reference), quote=F)

## pull out only ref-only and alt=SNP records
print(paste("Loading MIX",mixSampleName), quote=F)
mixSample = loadVCF(vcfTxtFile=paste("STAR-",reference,"-",mixSampleName,"/Aligned.sortedByCoord.f2.sorted.mpileup.call.txt",sep=""))
mixSample = mixSample[nchar(mixSample$ALT)<2,]

print(paste("Loading REF",refSampleName), quote=F)
refSample = loadVCF(vcfTxtFile=paste("STAR-",reference,"-",refSampleName,"/Aligned.sortedByCoord.f2.sorted.mpileup.call.txt",sep=""))
refSample = refSample[nchar(refSample$ALT)<2,]

print(paste("Loading ALT",altSampleName), quote=F)
altSample = loadVCF(vcfTxtFile=paste("STAR-",reference,"-",altSampleName,"/Aligned.sortedByCoord.f2.sorted.mpileup.call.txt",sep=""))
altSample = altSample[nchar(altSample$ALT)<2,]

## each data frame is
## CHR POS  REF ALT  REFfor REFrev  ALTfor ALTrev  REFtot ALTtot

merged = merge(altSample, refSample, by=c(1,2,3), suffixes=c(".a",".r"))
merged = merge(merged, mixSample, by=c(1,2,3,4))
merged$ALT.r = NULL

## CHR POS REF ALT.a REFfor.a REFrev.a ALTfor.a ALTrev.a REFtot.a ALTtot.a REFfor.r REFrev.r ALTfor.r ALTrev.r REFtot.r ALTtot.r REFfor REFrev ALTfor ALTrev REFtot ALTtot
colnames(merged) = c("CHR","POS","REF","ALT",
                     "REFfor.a","REFrev.a","ALTfor.a","ALTrev.a","REFtot.a","ALTtot.a",
                     "REFfor.r","REFrev.r","ALTfor.r","ALTrev.r","REFtot.r","ALTtot.r",
                     "REFfor.m","REFrev.m","ALTfor.m","ALTrev.m","REFtot.m","ALTtot.m")

rownames(merged) = paste(merged$CHR, merged$POS, sep=":")

## add in the total counts for each sample
merged$TOT.a = merged$REFtot.a + merged$ALTtot.a
merged$TOT.r = merged$REFtot.r + merged$ALTtot.r
merged$TOT.m = merged$REFtot.m + merged$ALTtot.m

## add in the ratio of mix/ref REF counts and mix/alt ALT counts
merged$REFratio = merged$REFtot.m/merged$REFtot.r
merged$ALTratio = merged$ALTtot.m/merged$ALTtot.a

## append the gene names (and update filters)
print("Appending gene names to loci", quote=F)
source("append-genes.R")

## create filters on final merged DF
print("Creating filters", quote=F)
source("filter-merged.R")
