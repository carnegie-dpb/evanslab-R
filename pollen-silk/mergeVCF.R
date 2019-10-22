## merge the data frames that are loaded with loadVCF
source("loadVCF.R")

reference = "B73"

refSampleName = "SS354+355"
altSampleName = "PS422"
mixSampleName = "S35464"

print(paste("Loading",refSampleName), quote=F)
refSample = loadVCF(vcfTxtFile=paste("STAR-",reference,"-",refSampleName,"/Aligned.sortedByCoord.f2.sorted.mpileup.call.txt",sep=""))
print(paste("Loading",altSampleName), quote=F)
altSample = loadVCF(vcfTxtFile=paste("STAR-",reference,"-",altSampleName,"/Aligned.sortedByCoord.f2.sorted.mpileup.call.txt",sep=""))
print(paste("Loading",mixSampleName), quote=F)
mixSample = loadVCF(vcfTxtFile=paste("STAR-",reference,"-",mixSampleName,"/Aligned.sortedByCoord.f2.sorted.mpileup.call.txt",sep=""))

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

## add in the total counts for each sample
merged$TOT.a = merged$REFtot.a + merged$ALTtot.a
merged$TOT.r = merged$REFtot.r + merged$ALTtot.r
merged$TOT.m = merged$REFtot.m + merged$ALTtot.m

## add in the ratio of mix/ref REF counts and mix/alt ALT counts
merged$REFratio = merged$REFtot.m/merged$REFtot.r
merged$ALTratio = merged$ALTtot.m/merged$ALTtot.a

## for "clean" version, keep rows which:
## 1. have nonzero alt sample ALT reads and ref sample REF reads
## 2. have alt sample REF reads<1% or ref sample ALT reads<1%
## 3. have nonzero mix sample REF reads
merged.clean = merged[
    merged$ALTtot.a>0 & merged$REFtot.r>0 &
    merged$REFtot.a/(merged$REFtot.a+merged$ALTtot.a)<0.01 & merged$ALTtot.r/(merged$REFtot.r+merged$ALTtot.r)<0.01 &
    merged$REFtot.m>0
   ,]

## for "strong" version, keep rows which:
## 1. have >10 total reads for all three samples
## 2. have >0 reads for ref sample REF reads and mix sample REF reads
merged.strong = merged[
    (merged$REFtot.a + merged$ALTtot.a)>0 &
    (merged$REFtot.r + merged$ALTtot.r)>0 &
    (merged$REFtot.m + merged$ALTtot.m)>0 &
    merged$REFtot.m>0 &
    merged$REFtot.r>0
   ,]

## for "bigref" version, keep rows which:
## 1. have >100 REF reads for ref sample
## 2. have >10  REF reads for mix sample
merged.bigref = merged[
    merged$REFtot.r>100 &
    merged$REFtot.m>10
   ,]


