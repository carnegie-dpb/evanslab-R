##
## load and merge the three B73-mapped files of interest
##

source("loadVCF.R")

## load
B73.YX24 = loadVCF("STAR-B73-YX24/Aligned.sortedByCoord.f2.sorted.mpileup.call.txt")
B73.S364_354 = loadVCF("STAR-B73-S364_354/Aligned.sortedByCoord.f2.sorted.mpileup.call.txt")
B73.SS364 = loadVCF("STAR-B73-SS364/Aligned.sortedByCoord.f2.sorted.mpileup.call.txt")

## merge
source("mergeVCF.R")
